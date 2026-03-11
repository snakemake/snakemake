from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Iterable, Tuple, Set

from sqlalchemy import (
    JSON,
    Boolean,
    Float,
    Integer,
    String,
    create_engine,
    select,
    delete,
    event,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, Session, mapped_column

from snakemake.persistence import MetadataRecord, PersistenceBase


class Base(DeclarativeBase):
    pass


class MetadataRecordORM(Base):
    __tablename__ = "snakemake_metadata"

    target: Mapped[str] = mapped_column(String, primary_key=True)

    incomplete: Mapped[bool] = mapped_column(Boolean, default=False)
    external_jobid: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    starttime: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    endtime: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    job_hash: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    rule: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    code: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    shellcmd: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    record_format_version: Mapped[int] = mapped_column(Integer, default=0)
    conda_env: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    container_img_url: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    software_stack_hash: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    input: Mapped[Optional[List[str]]] = mapped_column(JSON, nullable=True)
    log: Mapped[Optional[List[str]]] = mapped_column(JSON, nullable=True)
    params: Mapped[Optional[List[Any]]] = mapped_column(JSON, nullable=True)
    input_checksums: Mapped[Optional[Dict[str, Any]]] = mapped_column(
        JSON, default=dict
    )

    def to_record(self) -> MetadataRecord:
        return MetadataRecord(
            rule=self.rule,
            input=self.input,
            log=self.log,
            shellcmd=self.shellcmd,
            params=self.params,
            code=self.code,
            record_format_version=self.record_format_version,
            conda_env=self.conda_env,
            container_img_url=self.container_img_url,
            software_stack_hash=self.software_stack_hash,
            job_hash=self.job_hash,
            starttime=self.starttime,
            endtime=self.endtime,
            incomplete=self.incomplete,
            external_jobid=self.external_jobid,
            input_checksums=self.input_checksums or {},
        )


class LockORM(Base):
    __tablename__ = "snakemake_locks"
    file_path: Mapped[str] = mapped_column(String, primary_key=True)
    lock_type: Mapped[str] = mapped_column(String, primary_key=True)


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
        db_url: Optional[str] = None,
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

        self.engine = create_engine(db_url)

        @event.listens_for(self.engine, "connect")
        def set_sqlite_pragma(dbapi_connection, connection_record):
            if self.engine.dialect.name == "sqlite":
                cursor = dbapi_connection.cursor()
                # we may want to try this in the future if we encounter locking issues, but it can cause problems on network filesystems
                # cursor.execute("PRAGMA journal_mode=WAL")
                # cursor.execute("PRAGMA synchronous=NORMAL")
                cursor.execute("PRAGMA busy_timeout=10000")
                cursor.close()

        Base.metadata.create_all(self.engine)

    def _clear_cache(self) -> None:
        self._read_record_cached.cache_clear()

    @lru_cache(maxsize=None)
    def _read_record_cached(self, key: str) -> Optional[MetadataRecord]:
        with Session(self.engine) as session:
            record: Optional[MetadataRecordORM] = session.get(MetadataRecordORM, key)
            return record.to_record() if record else None

    def _read_record(self, key: str) -> Optional[MetadataRecord]:
        return self._read_record_cached(key)

    def _write_record(self, key: str, record: MetadataRecord) -> None:
        with Session(self.engine) as session:
            orm_record = session.get(MetadataRecordORM, key) or MetadataRecordORM(
                target=key
            )
            for field_name in record.keys():
                setattr(orm_record, field_name, getattr(record, field_name))
            session.add(orm_record)
            session.commit()

    def _delete_record(self, key: str) -> bool:
        with Session(self.engine) as session:
            if record := session.get(MetadataRecordORM, key):
                session.delete(record)
                session.commit()
                return True
        return False

    def _mark_incomplete(self, key: str, external_jobid: Optional[str]) -> None:
        with Session(self.engine) as session:
            record = session.get(MetadataRecordORM, key) or MetadataRecordORM(
                target=key
            )
            record.incomplete = True
            record.external_jobid = external_jobid
            session.add(record)
            session.commit()

    def _unmark_incomplete(self, key: str) -> None:
        with Session(self.engine) as session:
            if record := session.get(MetadataRecordORM, key):
                record.incomplete = False

                # if the record has no actual metadata (i.e., is just a stub from being marked incomplete),
                # delete it to correctly handle --drop-metadata and has_metadata() checks
                if record.record_format_version == 0:
                    session.delete(record)
                else:
                    session.add(record)

                session.commit()

    def _chunked(self, keys_list: List[str], size: int = 900) -> Iterable[List[str]]:
        for i in range(0, len(keys_list), size):
            yield keys_list[i : i + size]

    def _filter_incomplete_keys(self, keys: Iterable[str]) -> Set[str]:
        keys_list = list(keys)
        if not keys_list:
            return set()

        result = set()
        with Session(self.engine) as session:
            for chunk in self._chunked(keys_list):
                stmt = select(MetadataRecordORM.target).where(
                    MetadataRecordORM.target.in_(chunk),
                    MetadataRecordORM.incomplete == True,
                )
                result.update(session.scalars(stmt).all())
        return result

    def _get_external_jobids(self, keys: Iterable[str]) -> Set[str]:
        keys_list = list(keys)
        if not keys_list:
            return set()

        result = set()
        with Session(self.engine) as session:
            for chunk in self._chunked(keys_list):
                stmt = select(MetadataRecordORM.external_jobid).where(
                    MetadataRecordORM.target.in_(chunk),
                    MetadataRecordORM.external_jobid.is_not(None),
                )
                result.update(session.scalars(stmt).all())
        return result

    def _read_locks(self) -> Iterable[Tuple[str, str]]:
        with Session(self.engine) as session:
            stmt = select(LockORM)
            return [
                (lock.lock_type, lock.file_path) for lock in session.scalars(stmt).all()
            ]

    def _write_locks(self, lock_type: str, keys: Iterable[str]) -> None:
        keys_list = list(keys)
        if not keys_list:
            return

        with Session(self.engine) as session:
            for key in keys_list:
                session.add(LockORM(file_path=key, lock_type=lock_type))
            try:
                session.commit()
            except IntegrityError:
                session.rollback()
                raise snakemake.exceptions.LockException()

    def _delete_locks(self) -> None:
        with Session(self.engine) as session:
            session.execute(delete(LockORM))
            session.commit()
