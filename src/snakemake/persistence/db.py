from collections import OrderedDict
from pathlib import Path
from typing import Iterable

from sqlalchemy import (
    create_engine,
    select,
    delete,
    event,
)
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import DeclarativeBase, Session
from sqlmodel import SQLModel, Field

from snakemake.persistence import MetadataRecord, PersistenceBase
import snakemake.exceptions


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

        # ensure default db path is based on persistence path
        if db_url is None:
            db_url = f"sqlite:///{self.path / 'metadata.db'}"

        # use the absolute workdir path as a namespace
        # to allow using the same db for multiple different Snakemake instances running in different directories
        self.namespace = str(self.path.absolute())

        # using custom LRU to be able to clear specific keys only, see https://stackoverflow.com/a/52101715
        self._metadata_cache = OrderedDict()
        self._cache_size = 16384

        self.engine = create_engine(db_url)

        # for SQLite, set a busy timeout to avoid immediate failures on database locks;
        # if available, use the latency_wait setting instead.
        # TODO: if that doesn't help, consider using flufl.lock
        #  or a NullPool or enabling BEGIN IMMEDIATE
        busy_timeout = 10000
        if self.dag:
            latency_wait_s = self.dag.workflow.execution_settings.latency_wait * 1000
            busy_timeout = max(busy_timeout, latency_wait_s)

        @event.listens_for(self.engine, "connect")
        def set_sqlite_pragma(dbapi_connection, connection_record):
            if self.engine.dialect.name == "sqlite":
                cursor = dbapi_connection.cursor()
                # we may want to try this in the future if we encounter locking issues, but it can cause problems on network filesystems
                # cursor.execute("PRAGMA journal_mode=WAL")
                # cursor.execute("PRAGMA synchronous=NORMAL")
                cursor.execute(f"PRAGMA busy_timeout={busy_timeout}")
                cursor.close()

        SQLModel.metadata.create_all(self.engine)

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
