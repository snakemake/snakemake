import asyncio
import os
import time
from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

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

from snakemake.persistence import MetadataRecord, PersistenceBase, RECORD_FORMAT_VERSION
from snakemake_interface_executor_plugins.settings import ExecMode
import snakemake.exceptions
from snakemake.logging import logger
from snakemake.io import IOCache, _IOFile


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
    lock_type: Mapped[str] = mapped_column(String)


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
        db_url: str = "sqlite:///.snakemake/metadata.db",
    ):
        if path is None:
            self._path = Path(os.path.abspath(".snakemake"))
        else:
            self._path = path
        os.makedirs(self.path, exist_ok=True)

        self.dag = dag

        self.conda_env_archive_path = os.path.join(self.path, "conda-archive")
        self.benchmark_path = os.path.join(self.path, "benchmarks")
        self.source_cache = os.path.join(self.path, "source_cache")
        self.iocache_path = os.path.join(self.path, "iocache")
        # place to store any auxiliary information needed during a run (e.g. source tarballs)
        self._aux_path = os.path.join(self.path, "auxiliary")

        if conda_prefix is None:
            self.conda_env_path = os.path.join(self.path, "conda")
        else:
            self.conda_env_path = os.path.abspath(conda_prefix)
        if singularity_prefix is None:
            self.container_img_path = os.path.join(self.path, "singularity")
        else:
            self.container_img_path = os.path.abspath(singularity_prefix)
        if shadow_prefix is None:
            self.shadow_path = os.path.join(self.path, "shadow")
        else:
            self.shadow_path = os.path.join(shadow_prefix, "shadow")

        for d in (
            self.shadow_path,
            self.conda_env_archive_path,
            self.conda_env_path,
            self.container_img_path,
            self.aux_path,
            self.iocache_path,
        ):
            os.makedirs(d, exist_ok=True)

        if nolock:
            self.lock = self.noop
            self.unlock = self.noop
        if warn_only:
            self.lock = self.lock_warn_only
            self.unlock = self.noop

        self.max_checksum_file_size = (
            self.dag.workflow.dag_settings.max_checksum_file_size
        )

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

    @property
    def path(self) -> Path:
        return Path(self._path)

    @property
    def aux_path(self) -> Path:
        return Path(self._aux_path)

    def _get_key(self, f: _IOFile) -> str:
        assert isinstance(f, _IOFile)
        return str(f.storage_object.query if f.is_storage else f)

    @lru_cache(maxsize=None)
    def _read_record_cached(self, key: str) -> Optional[MetadataRecord]:
        with Session(self.engine) as session:
            record: Optional[MetadataRecordORM] = session.get(MetadataRecordORM, key)
            if record:
                return record.to_record()
        return None

    def metadata(self, target: Any) -> Optional[MetadataRecord]:
        return self._read_record_cached(self._get_key(target))

    def cleanup_metadata(self, target: Any) -> bool:
        self._read_record_cached.cache_clear()

        key = self._get_key(target)
        with Session(self.engine) as session:
            record = session.get(MetadataRecordORM, key)
            if record:
                session.delete(record)
                session.commit()
                return True
        return False

    def started(self, job: Any, external_jobid: Optional[str] = None) -> None:
        self._read_record_cached.cache_clear()
        with Session(self.engine) as session:
            for target in job.output:
                key = self._get_key(target)
                record = session.get(MetadataRecordORM, key) or MetadataRecordORM(
                    target=key
                )
                record.incomplete = True
                record.external_jobid = external_jobid
                session.add(record)
            session.commit()

    async def finished(self, job: Any) -> None:
        self._read_record_cached.cache_clear()

        if not self.dag.workflow.execution_settings.keep_metadata:
            with Session(self.engine) as session:
                for f in job.output:
                    record = session.get(MetadataRecordORM, self._get_key(f))
                    if record:
                        record.incomplete = False
                        session.add(record)
                session.commit()
            return

        if (
            self.dag.workflow.exec_mode == ExecMode.DEFAULT
            or self.dag.workflow.remote_execution_settings.immediate_submit
        ):
            code = self._code(job.rule)
            input_files = self._input(job)
            log = self._log(job)
            params = self._params(job)
            shellcmd = job.shellcmd
            conda_env = self._conda_env(job)
            software_stack_hash = self._software_stack_hash(job)
            fallback_time = time.time()
            job_hash_val = hash(job)

            with Session(self.engine) as session:
                for f in job.output:
                    key = self._get_key(f)

                    endtime = (
                        (await f.mtime()).local_or_storage()
                        if await f.exists()
                        else fallback_time
                    )
                    checksums = {
                        str(infile): await infile.checksum(self.max_checksum_file_size)
                        for infile in job.input
                    }
                    checksums = {k: v for k, v in checksums.items() if v is not None}

                    record = session.get(MetadataRecordORM, key) or MetadataRecordORM(
                        target=key
                    )

                    if record.starttime is None:
                        record.starttime = fallback_time

                    record.incomplete = False
                    record.record_format_version = RECORD_FORMAT_VERSION
                    record.rule = job.rule.name
                    record.code = code
                    record.input = input_files
                    record.log = log
                    record.params = params
                    record.shellcmd = shellcmd
                    record.endtime = endtime
                    record.job_hash = job_hash_val
                    record.conda_env = conda_env
                    record.software_stack_hash = software_stack_hash
                    record.container_img_url = job.container_img_url
                    record.input_checksums = checksums

                    session.add(record)
                session.commit()

    async def incomplete(self, job: Any) -> List[Any]:
        output_ids = [self._get_key(f) for f in job.output]

        with Session(self.engine) as session:
            stmt = select(MetadataRecordORM.target).where(
                MetadataRecordORM.target.in_(output_ids),
                MetadataRecordORM.incomplete == True,
            )
            incomplete_paths_in_db = set(session.scalars(stmt).all())

        async def is_incomplete(f):
            exists = await f.exists()
            marked = self._get_key(f) in incomplete_paths_in_db
            return f if exists and marked else None

        async with asyncio.TaskGroup() as tg:
            tasks = [tg.create_task(is_incomplete(f)) for f in job.output]

        return [task.result() for task in tasks if task.result() is not None]

    def external_jobids(self, job: Any) -> List[str]:
        output_ids = [self._get_key(f) for f in job.output]
        with Session(self.engine) as session:
            stmt = select(MetadataRecordORM.external_jobid).where(
                MetadataRecordORM.target.in_(output_ids),
                MetadataRecordORM.external_jobid.is_not(None),
            )
            return list(set(session.scalars(stmt).all()))

    @property
    def locked(self) -> bool:
        inputfiles = set(str(f) for f in self.all_inputfiles())
        outputfiles = set(str(f) for f in self.all_outputfiles())

        with Session(self.engine) as session:
            stmt = select(LockORM)
            locks = session.scalars(stmt).all()
            for lock in locks:
                if lock.lock_type == "input" and lock.file_path in outputfiles:
                    return True
                if lock.lock_type == "output" and (
                    lock.file_path in outputfiles or lock.file_path in inputfiles
                ):
                    return True
        return False

    @contextmanager
    def lock_warn_only(self):
        if self.locked:
            logger.info(
                "Error: Directory cannot be locked. Another Snakemake instance is running."
            )
        yield

    @contextmanager
    def lock(self):
        if self.locked:
            raise snakemake.exceptions.LockException()
        try:
            with Session(self.engine) as session:
                for f in self.all_inputfiles():
                    session.merge(LockORM(file_path=str(f), lock_type="input"))
                for f in self.all_outputfiles():
                    session.merge(LockORM(file_path=str(f), lock_type="output"))
                session.commit()
            yield
        finally:
            self.unlock()

    def unlock(self) -> None:
        logger.debug("unlocking database locks")
        with Session(self.engine) as session:
            session.execute(delete(LockORM))
            session.commit()

    def cleanup_locks(self) -> None:
        self.unlock()

    @contextmanager
    def noop(self, *args):
        yield

    def deactivate_cache(self) -> None:
        self._read_record_cached.cache_clear()

    @property
    def _iocache_filename(self):
        return os.path.join(self.iocache_path, "latest.pkl")

    def save_iocache(self) -> None:
        with open(self._iocache_filename, "wb") as handle:
            self.dag.workflow.iocache.save(handle)

    def load_iocache(self) -> None:
        if os.path.exists(self._iocache_filename):
            logger.info("Loading trusted IOCache from latest dry-run.")
            with open(self._iocache_filename, "rb") as handle:
                self.dag.workflow.iocache = IOCache.load(handle)

    def drop_iocache(self) -> None:
        if os.path.exists(self._iocache_filename):
            os.remove(self._iocache_filename)
