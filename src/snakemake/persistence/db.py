import functools
import os
from collections import OrderedDict
from pathlib import Path
from typing import Iterable

import psutil
from sqlalchemy import (
    create_engine,
    select,
    delete,
    event,
)
from sqlalchemy.engine.url import make_url, URL
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import DeclarativeBase, Session
from sqlmodel import SQLModel, Field
import sqlite3

from snakemake.persistence import MetadataRecord, PersistenceBase
from snakemake.logging import logger
import snakemake.exceptions


class SettingsORM(SQLModel, table=True):
    __tablename__ = "snakemake_persistence_db_settings"
    namespace: str = Field(primary_key=True)
    is_network_fs: bool


class Base(DeclarativeBase):
    pass


class MetadataRecordORM(MetadataRecord, table=True):
    __tablename__ = "snakemake_metadata"
    namespace: str = Field(primary_key=True)
    target: str = Field(primary_key=True)


class LockORM(SQLModel, table=True):
    __tablename__ = "snakemake_locks"
    namespace: str = Field(primary_key=True)
    file_path: str = Field(primary_key=True)
    lock_type: str = Field(primary_key=True)


class DbPersistence(PersistenceBase):
    def __init__(
        self,
        nolock=False,
        dag=None,
        conda_prefix=None,
        singularity_prefix=None,
        shadow_prefix=None,
        warn_only=False,
        path: Path | None = None,
        db_url: str | None = None,
    ):
        super().__init__(
            nolock=nolock,
            dag=dag,
            conda_prefix=conda_prefix,
            singularity_prefix=singularity_prefix,
            shadow_prefix=shadow_prefix,
            warn_only=warn_only,
            path=path,
        )

        # use the absolute workdir path as a namespace
        # to allow using the same db for multiple different Snakemake instances running in different directories
        self.namespace = str(self.path.absolute())

        # using custom LRU to be able to clear specific keys only, see https://stackoverflow.com/a/52101715
        self._metadata_cache = OrderedDict()
        self._cache_size = 16384

        self._setup_database(db_url)

    def _setup_database(self, db_url: str | None) -> None:
        """Sets up the database connection and schema."""

        # ensure default db path is based on persistence path
        if db_url is None:
            db_url = f"sqlite:///{self.path / 'metadata.db'}"

        parsed_url = make_url(db_url)
        is_network_fs = self._check_network_fs(parsed_url)
        busy_timeout = self._get_busy_timeout()

        self.engine = create_engine(db_url)
        self._attach_sqlite_pragmas(is_network_fs, busy_timeout)

        SQLModel.metadata.create_all(self.engine)

    def _check_network_fs(self, parsed_url: URL) -> bool:
        """
        Determines if the sqlite DB is on a network FS, caching the result in the DB.
        Tries looking up the result in the snakemake_persistence_db_settings table,
        otherwise executes psutil to figure out the result.
        Always returns False for non-SQLite URLs.
        """
        if parsed_url.get_backend_name() != "sqlite":
            return False

        sqlite_db_path = parsed_url.database
        if not sqlite_db_path or sqlite_db_path == ":memory:":
            return False

        # raw connection to avoid PRAGMA chicken-vs-egg problems
        with sqlite3.connect(sqlite_db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(
                "CREATE TABLE IF NOT EXISTS snakemake_persistence_db_settings "
                "(namespace TEXT PRIMARY KEY, is_network_fs BOOLEAN)"
            )
            cursor.execute(
                "SELECT is_network_fs FROM snakemake_persistence_db_settings WHERE namespace = ?",
                (self.namespace,),
            )
            row = cursor.fetchone()

            if row is not None:
                return bool(row[0])

            is_network_fs = is_network_filesystem(sqlite_db_path)
            cursor.execute(
                "INSERT OR IGNORE INTO snakemake_persistence_db_settings (namespace, is_network_fs) VALUES (?, ?)",
                (self.namespace, is_network_fs),
            )
            conn.commit()

            return is_network_fs

    def _get_busy_timeout(self) -> int:
        """
        Determines SQLite busy timeout (`max(10000, latency_wait * 1000)`).
        """
        base_timeout = 10000
        if self.dag:
            latency_wait_s = self.dag.workflow.execution_settings.latency_wait * 1000
            return max(base_timeout, latency_wait_s)
        return base_timeout

    def _attach_sqlite_pragmas(self, is_network_fs: bool, busy_timeout: int) -> None:
        """Binds the PRAGMA execution to the SQLAlchemy connection event."""
        if self.engine.dialect.name != "sqlite":
            return

        @event.listens_for(self.engine, "connect")
        def set_sqlite_pragma(dbapi_connection, connection_record):
            cursor = dbapi_connection.cursor()

            if is_network_fs:
                cursor.execute("PRAGMA journal_mode=PERSIST")
                cursor.execute("PRAGMA synchronous=OFF")
                cursor.execute("PRAGMA temp_store=MEMORY")
                cursor.execute("PRAGMA cache_size=-64000")
            else:
                cursor.execute("PRAGMA journal_mode=TRUNCATE")
                cursor.execute("PRAGMA synchronous=NORMAL")

            cursor.execute(f"PRAGMA busy_timeout={busy_timeout}")
            cursor.close()

    def _clear_cache(self) -> None:
        self._metadata_cache.clear()

    def _invalidate_cache(self, key: str) -> None:
        self._metadata_cache.pop(key, None)

    def _read_record(self, key: str) -> MetadataRecord | None:
        if key in self._metadata_cache:
            self._metadata_cache.move_to_end(key)
            return self._metadata_cache[key]

        with Session(self.engine) as session:
            record: MetadataRecordORM | None = session.get(
                MetadataRecordORM, (self.namespace, key)
            )
            record_ = MetadataRecord.model_validate(record) if record else None

            self._metadata_cache[key] = record_
            if len(self._metadata_cache) > self._cache_size:
                self._metadata_cache.popitem(last=False)

            return record_

    def _write_record(self, key: str, record: MetadataRecord) -> None:
        self._invalidate_cache(key)
        with Session(self.engine) as session:
            orm_record = session.get(MetadataRecordORM, (self.namespace, key))
            if orm_record:
                orm_record.sqlmodel_update(record)
            else:
                orm_record = MetadataRecordORM(
                    namespace=self.namespace, target=key, **record.model_dump()
                )

            session.add(orm_record)
            session.commit()

    def _delete_record(self, key: str) -> bool:
        self._invalidate_cache(key)
        with Session(self.engine) as session:
            if record := session.get(MetadataRecordORM, (self.namespace, key)):
                session.delete(record)
                session.commit()
                return True
        return False

    def _mark_incomplete(self, key: str, external_jobid: str | None) -> None:
        self._invalidate_cache(key)
        with Session(self.engine) as session:
            record = session.get(
                MetadataRecordORM, (self.namespace, key)
            ) or MetadataRecordORM(namespace=self.namespace, target=key)
            record.incomplete = True
            record.external_jobid = external_jobid
            session.add(record)
            session.commit()

    def _unmark_incomplete(self, key: str) -> None:
        self._invalidate_cache(key)
        with Session(self.engine) as session:
            if record := session.get(MetadataRecordORM, (self.namespace, key)):
                record.incomplete = False

                # if the record has no actual metadata (i.e., is just a stub from being marked incomplete),
                # delete it to correctly handle --drop-metadata and has_metadata() checks
                if record.record_format_version == 0:
                    session.delete(record)
                else:
                    session.add(record)

                session.commit()

    def _chunked(self, keys_list: list[str], size: int = 900) -> Iterable[list[str]]:
        for i in range(0, len(keys_list), size):
            yield keys_list[i : i + size]

    def _filter_incomplete_keys(self, keys: Iterable[str]) -> set[str]:
        if not (keys_list := list(keys)):
            return set()

        result = set()
        with Session(self.engine) as session:
            for chunk in self._chunked(keys_list):
                stmt = select(MetadataRecordORM.target).where(
                    MetadataRecordORM.namespace == self.namespace,
                    MetadataRecordORM.target.in_(chunk),
                    MetadataRecordORM.incomplete == True,
                )
                result.update(session.scalars(stmt).all())
        return result

    def _get_external_jobids(self, keys: Iterable[str]) -> set[str]:
        if not (keys_list := list(keys)):
            return set()

        result = set()
        with Session(self.engine) as session:
            for chunk in self._chunked(keys_list):
                stmt = select(MetadataRecordORM.external_jobid).where(
                    MetadataRecordORM.namespace == self.namespace,
                    MetadataRecordORM.target.in_(chunk),
                    MetadataRecordORM.external_jobid.is_not(None),
                )
                result.update(session.scalars(stmt).all())
        return result

    def _read_locks(self) -> Iterable[tuple[str, str]]:
        with Session(self.engine) as session:
            stmt = select(LockORM).where(LockORM.namespace == self.namespace)
            return [
                (lock.lock_type, lock.file_path) for lock in session.scalars(stmt).all()
            ]

    def _write_locks(self, lock_type: str, keys: Iterable[str]) -> None:
        if not (keys_list := list(set(keys))):
            return

        with Session(self.engine) as session:
            for key in keys_list:
                session.add(
                    LockORM(
                        namespace=self.namespace, file_path=key, lock_type=lock_type
                    )
                )
            try:
                session.commit()
            except IntegrityError as e:
                session.rollback()
                raise snakemake.exceptions.LockException() from e

    def _delete_locks(self) -> None:
        with Session(self.engine) as session:
            session.execute(delete(LockORM).where(LockORM.namespace == self.namespace))
            session.commit()


@functools.cache
def is_network_filesystem(path: Path | str) -> bool:
    """
    Detects if a given path resides on a network filesystem.
    """
    if os.name != "posix":
        return False

    path_obj = Path(path).resolve()

    logger.info("Determining filesystem type for metadata persistence storage...")
    try:
        best_match = max(
            (
                mount
                for mount in psutil.disk_partitions(all=True)
                if path_obj.is_relative_to(Path(mount.mountpoint).resolve())
            ),
            key=lambda part: len(Path(part.mountpoint).resolve().parts),
            default=None,
        )
        fstype = best_match.fstype if best_match else ""

        network_fs_types = {
            "afs",
            "beegfs",
            "ceph",
            "cifs",
            "fhgfs",
            "fuse.juicefs",
            "fuse.sshfs",
            "glusterfs",
            "gpfs",
            "lustre",
            "nfs",
            "nfs3",
            "nfs4",
            "orangefs",
            "pvfs2",
            "smbfs",
        }
        return fstype.casefold() in network_fs_types

    except (PermissionError, OSError):
        return False
