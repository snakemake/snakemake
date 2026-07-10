from snakemake.common import is_serializable
import asyncio
import os
import shutil
import time
from abc import abstractmethod
from contextlib import contextmanager
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import (
    Any,
    Callable,
    Iterable,
)

from snakemake_interface_executor_plugins.persistence import (
    PersistenceExecutorInterface,
)
from snakemake_interface_executor_plugins.settings import ExecMode
from sqlalchemy import Column, JSON
from sqlmodel import SQLModel, Field

from snakemake.io import get_flag_value, is_flagged, _IOFile, IOCache
from snakemake.logging import logger
import snakemake.exceptions

RECORD_FORMAT_VERSION = 7
UNREPRESENTABLE = object()


class MetadataRecord(SQLModel):
    rule: str | None = None
    input: list[str] | None = Field(default=None, sa_column=Column(JSON))
    log: list[str] | None = Field(default=None, sa_column=Column(JSON))
    shellcmd: str | None = None
    params: list[Any] | None = Field(default=None, sa_column=Column(JSON))
    code: str | None = None
    record_format_version: int = 0
    software_stack_hash: str | None = None
    software: list[str] | None = Field(default=None, sa_column=Column(JSON))
    job_hash: int | None = None
    starttime: float | None = None
    endtime: float | None = None
    incomplete: bool | None = None
    external_jobid: str | None = None
    input_checksums: dict[str, Any] | None = Field(
        default_factory=dict, sa_column=Column(JSON)
    )

    def __getitem__(self, key: str) -> Any:
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(key)

    def get(self, key: str, default: Any = None) -> Any:
        return getattr(self, key, default)

    def keys(self):
        return self.model_dump().keys()

    def items(self):
        return self.model_dump().items()


@dataclass
class ParamsChange:
    only_old: set[Any] = field(default_factory=set)
    only_new: set[Any] = field(default_factory=set)
    files: set[str] = field(default_factory=set)

    def __bool__(self):
        return bool(self.only_old or self.only_new)

    def __or__(self, other):
        if not self:
            return other
        if not other:
            return self
        return ParamsChange(
            only_old=self.only_old | other.only_old,
            only_new=self.only_new | other.only_new,
            files=self.files | other.files,
        )

    def __iter__(self):
        return iter(self.files)

    def __str__(self):
        if not self:
            return "No params change"
        else:

            def fmt_set(s, label):
                if s:
                    return f"{label}: {','.join(s)}"
                else:
                    return f"{label}: <nothing exclusive>"

            return (
                "Union of exclusive params before and now across all output: "
                f"{fmt_set(self.only_old, 'before')} "
                f"{fmt_set(self.only_new, 'now')} "
            )


NO_PARAMS_CHANGE = ParamsChange()


def _bool_or_gen(func, job, file=None):
    if file is None:
        return (f for f in job.output if func(job, file=f))
    else:
        return func(job, file=file)


class EnvironmentMaintenanceMixin:
    """Handles standard filesystem cleanups (shadow, software deployment)."""

    shadow_path: str | Path
    container_img_path: str | Path
    conda_env_archive_path: str | Path
    dag: Any

    def cleanup_shadow(self) -> None:
        if os.path.exists(self.shadow_path):
            shutil.rmtree(self.shadow_path)
            os.mkdir(self.shadow_path)


class FileBackupMixin:
    """Handles backing up and restoring output files locally."""

    path: Path

    @property
    def _backup_path(self) -> Path:
        return self.path / "backups"

    def _get_backup_path(self, path: Path) -> Path:
        if path.is_absolute():
            path = path.relative_to(path.parents[-1])
        return self._backup_path / path

    def backup_output(self, path: Path) -> None:
        backup_path = self._get_backup_path(path)
        backup_path.parent.mkdir(parents=True, exist_ok=True)
        if path.is_dir():
            shutil.copytree(path, backup_path)
        else:
            shutil.copy(path, backup_path)

    def restore_output(self, path: Path) -> bool:
        backup_path = self._get_backup_path(path)
        if not backup_path.exists():
            logger.warning(f"Cannot restore {path}: no backup found.")
            return False
        elif backup_path.is_dir():
            if path.exists():
                shutil.rmtree(path)
            shutil.copytree(backup_path, path)
            shutil.rmtree(backup_path)
            return True
        else:
            shutil.copy(backup_path, path)
            backup_path.unlink()
            return True

    def cleanup_backup(self, path: Path) -> None:
        backup_path = self._get_backup_path(path)
        if backup_path.exists():
            if backup_path.is_dir():
                shutil.rmtree(backup_path)
            else:
                backup_path.unlink()


class PersistenceBase(
    PersistenceExecutorInterface, EnvironmentMaintenanceMixin, FileBackupMixin
):
    dag: Any
    max_checksum_file_size: int

    def __init__(
        self,
        nolock=False,
        dag=None,
        shadow_prefix=None,
        warn_only=False,
        path: Path | None = None,
    ):
        if path is None:
            self._path = Path(os.path.abspath(".snakemake"))
        else:
            self._path = path
        os.makedirs(self.path, exist_ok=True)

        self.dag = dag

        self.benchmark_path = os.path.join(self.path, "benchmarks")
        self.source_cache = os.path.join(self.path, "source_cache")
        self.iocache_path = os.path.join(self.path, "iocache")
        self._aux_path = os.path.join(self.path, "auxiliary")

        if shadow_prefix is None:
            self.shadow_path = os.path.join(self.path, "shadow")
        else:
            self.shadow_path = os.path.join(shadow_prefix, "shadow")

        for d in (
            self.shadow_path,
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

        if self.dag is not None:
            self.max_checksum_file_size = (
                self.dag.workflow.dag_settings.max_checksum_file_size
            )
        else:
            self.max_checksum_file_size = 0

    @property
    def path(self) -> Path:
        return Path(self._path)

    @property
    def aux_path(self) -> Path:
        return Path(self._aux_path)

    @abstractmethod
    def _read_record(self, key: str) -> MetadataRecord | None:
        """
        Fetch the metadata record for the given key.

        Args:
            key (str): The unique string identifier for the output file.

        Returns:
            MetadataRecord | None: The deserialized record, or None if it does not exist.
        """
        ...

    @abstractmethod
    def _write_record(self, key: str, record: MetadataRecord) -> None:
        """
        Write or overwrite the metadata record for the given key.

        Args:
            key (str): The unique string identifier for the output file.
            record (MetadataRecord): The complete metadata record to store.
        """
        ...

    @abstractmethod
    def _delete_record(self, key: str) -> bool:
        """
        Delete the metadata record associated with the given key.

        Args:
            key (str): The unique string identifier for the output file.

        Returns:
            bool: True if the record was found and deleted, False if it did not exist.
        """
        ...

    @abstractmethod
    def _mark_incomplete(self, key: str, external_jobid: str | None) -> None:
        """
        Mark a file key as currently incomplete (job has started but not finished).

        Args:
            key (str): The unique string identifier for the output file.
            external_jobid (str | None): The cluster/executor job ID associated with this run.
        """
        ...

    @abstractmethod
    def _unmark_incomplete(self, key: str) -> None:
        """
        Remove the incomplete marker for a given file key.

        Args:
            key (str): The unique string identifier for the output file.
        """
        ...

    @abstractmethod
    def _filter_incomplete_keys(self, keys: Iterable[str]) -> set[str]:
        """
        Given a list of file keys, return a subset of those keys that are marked as incomplete.

        Args:
            keys (Iterable[str]): A collection of file keys to check.

        Returns:
            set[str]: The keys from the input iterable that have an active incomplete marker.
        """
        ...

    @abstractmethod
    def _get_external_jobids(self, keys: Iterable[str]) -> set[str]:
        """
        Given a list of file keys, return all unique external job IDs associated with them.

        Args:
            keys (Iterable[str]): A collection of file keys to check.

        Returns:
            set[str]: A set of all non-null external job IDs found for the provided keys.
        """
        ...

    @abstractmethod
    def _read_locks(self) -> Iterable[tuple[str, str]]:
        """
        Retrieve all currently active locks.

        Returns:
            Iterable[tuple[str, str]]: An iterable yielding (lock_type, key) pairs,
                                       where lock_type is 'input' or 'output'.
        """
        ...

    @abstractmethod
    def _write_locks(self, lock_type: str, keys: Iterable[str]) -> None:
        """
        Create locks of a specific type for multiple keys.

        Args:
            lock_type (str): The type of lock ('input' or 'output').
            keys (Iterable[str]): The file keys to lock.
        """
        ...

    @abstractmethod
    def _delete_locks(self) -> None:
        """
        Remove all active locks managed by this persistence instance.
        """
        ...

    @abstractmethod
    def _clear_cache(self) -> None:
        """
        Clear any internal, memory-based caches used by the persistence backend.
        This is called prior to state mutations to ensure data consistency.
        """
        ...

    def _get_key(self, f: _IOFile) -> str:
        assert isinstance(f, _IOFile)
        return str(f.storage_object.query if f.is_storage else f)

    def metadata(self, target: Any) -> MetadataRecord | None:
        return self._read_record(self._get_key(target))

    def cleanup_metadata(self, target: Any) -> bool:
        self._clear_cache()
        key = self._get_key(target)
        self._unmark_incomplete(key)
        return self._delete_record(key)

    def started(self, job: Any, external_jobid: str | None = None) -> None:
        self._clear_cache()
        for f in job.output:
            self._mark_incomplete(self._get_key(f), external_jobid)

    async def finished(self, job: Any) -> None:
        self._clear_cache()

        if not self.dag.workflow.execution_settings.keep_metadata:
            for f in job.output:
                self._unmark_incomplete(self._get_key(f))
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
            fallback_time = time.time()
            job_hash_val = hash(job)

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

                record = self._read_record(key) or MetadataRecord()

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
                record.input_checksums = checksums
                record.software = (
                    [
                        f"{software.name} {software.version}"
                        for software in job.software_env.report_software()
                        if not software.is_secondary
                    ]
                    if job.software_env is not None
                    else None
                )
                record.software_stack_hash = (
                    job.software_env.hash() if job.software_env else None
                )

                self._write_record(key, record)

        for f in job.output:
            self._unmark_incomplete(self._get_key(f))

    async def incomplete(self, job: Any) -> list[Any]:
        keys = [self._get_key(f) for f in job.output]
        marked_incomplete = self._filter_incomplete_keys(keys)

        if not marked_incomplete:
            return []

        async def is_incomplete(f):
            if self._get_key(f) in marked_incomplete and await f.exists():
                return f
            return None

        async with asyncio.TaskGroup() as tg:
            tasks = [tg.create_task(is_incomplete(f)) for f in job.output]

        return [t.result() for t in tasks if t.result() is not None]

    def external_jobids(self, job: Any) -> list[str]:
        keys = [self._get_key(f) for f in job.output]
        return list(self._get_external_jobids(keys))

    @property
    def locked(self) -> bool:
        inputfiles = set(str(f) for f in self.all_inputfiles())
        outputfiles = set(str(f) for f in self.all_outputfiles())

        for lock_type, key in self._read_locks():
            if lock_type == "input" and key in outputfiles:
                return True
            if lock_type == "output" and (key in outputfiles or key in inputfiles):
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
            self._write_locks("input", (str(f) for f in self.all_inputfiles()))
            self._write_locks("output", (str(f) for f in self.all_outputfiles()))
            yield
        finally:
            self.unlock()

    def unlock(self) -> None:
        logger.debug("unlocking")
        self._delete_locks()

    def cleanup_locks(self) -> None:
        self.unlock()

    def deactivate_cache(self) -> None:
        self._clear_cache()

    def cleanup(self, job: Any) -> None:
        for f in job.output:
            self.cleanup_metadata(f)

    def has_metadata(self, job: Any) -> bool:
        return all(self.metadata(path) is not None for path in job.output)

    def has_outdated_metadata(self, job: Any) -> bool:
        return any(
            (m := self.metadata(path)) is not None
            and m.record_format_version < RECORD_FORMAT_VERSION
            for path in job.output
        )

    def rule(self, path: Any) -> str | None:
        return (m := self.metadata(path)) and m.rule

    def input(self, path: Any) -> list[str] | None:
        return (m := self.metadata(path)) and m.input

    def log(self, path: Any) -> list[str] | None:
        return (m := self.metadata(path)) and m.log

    def shellcmd(self, path: Any) -> str | None:
        return (m := self.metadata(path)) and m.shellcmd

    def params(self, path: Any) -> list[Any] | None:
        return (m := self.metadata(path)) and m.params

    def code(self, path: Any) -> str | None:
        return (m := self.metadata(path)) and m.code

    def record_format_version(self, path: Any) -> int | None:
        return (m := self.metadata(path)) and m.record_format_version

    def conda_env(self, path: Any) -> str | None:
        return (m := self.metadata(path)) and m.conda_env

    def container_img_url(self, path: Any) -> str | None:
        return (m := self.metadata(path)) and m.container_img_url

    def software_stack_hash(self, path: Any) -> str | None:
        return (m := self.metadata(path)) and m.software_stack_hash

    def input_checksums(self, job: Any, input_path: Any) -> set[Any]:
        def ensure_algorithm(checksum):
            return (
                checksum
                if checksum is None or ":" in checksum
                else f"sha256:{checksum}"
            )

        return set(
            (
                ensure_algorithm(m.input_checksums.get(input_path))
                if (m := self.metadata(output_path)) and m.input_checksums
                else None
            )
            for output_path in job.output
        )

    def code_changed(self, job: Any, file=None) -> Any:
        return _bool_or_gen(self._code_changed, job, file=file)

    def input_changed(self, job: Any, file=None) -> Any:
        return _bool_or_gen(self._input_changed, job, file=file)

    def software_stack_changed(self, job: Any, file=None) -> Any:
        return _bool_or_gen(self._software_stack_changed, job, file=file)

    def _code_changed(self, job, file=None) -> bool:
        fmt_version = self.record_format_version(file)
        if fmt_version is None or fmt_version < 3:
            return False
        recorded = self.code(file)
        return recorded is not None and recorded != self._code(job.rule)

    def _input_changed(self, job, file=None) -> bool:
        fmt_version = self.record_format_version(file)
        if fmt_version is None or fmt_version < 4:
            return False
        recorded = self.input(file)
        return recorded is not None and recorded != self._input(job)

    def _software_stack_changed(self, job, file=None) -> bool:
        fmt_version = self.record_format_version(file)
        if fmt_version is None or fmt_version < 5:
            # no reliable software stack hash stored (previous storage ignored pin files
            # and aux deploy files of conda envs as well as env modules)
            return False
        recorded = self.software_stack_hash(file)
        current = job.software_env.hash() if job.software_env is not None else None
        return recorded is not None and recorded != current

    def params_changed(self, job: Any, file=None) -> bool | ParamsChange:
        files = [file] if file is not None else job.output
        changes = NO_PARAMS_CHANGE
        new = set(self._params(job))

        for outfile in files:
            fmt_version = self.record_format_version(outfile)
            if fmt_version is None or fmt_version < 6:
                # no reliable params stored (version 4 refactored params storage
                # and version 6 fixed a bug in determination of whether params are
                # derived from e.g. input or output files). If they are,
                # there is the risk to store storage paths here. Derived param
                # changes will also be captured by input changes.
                continue
            recorded = self.params(outfile)
            if recorded is not None:
                old = set(recorded)
                changes |= ParamsChange(
                    only_old=old - new, only_new=new - old, files={outfile}
                )
        return changes

    @lru_cache()
    def _code(self, rule) -> str | None:
        # Scripts and notebooks are triggered by changes in the script mtime.
        # Changes to python and shell rules are triggered by changes in the plain text.
        if rule.shellcmd is not None:
            return rule.shellcmd
        if rule.run_func_src is not None:
            return rule.run_func_src
        return None

    @lru_cache()
    def _input(self, job) -> list[str]:
        def get_paths():
            for f in job.input:
                if f.is_storage:
                    yield f.storage_object.query
                elif is_flagged(f, "pipe"):
                    yield "<pipe>"
                elif is_flagged(f, "service"):
                    yield "<service>"
                else:
                    yield (
                        # get the true path instead of the cache path
                        get_flag_value(f, "sourcecache_entry")
                        if is_flagged(f, "sourcecache_entry")
                        else f
                    )

        return sorted(get_paths())

    @lru_cache()
    def _log(self, job) -> list[str]:
        return sorted(job.log)

    def _serialize_param_builtin(self, value: Any) -> str | object:
        if is_serializable(value):
            return repr(value)
        else:
            return UNREPRESENTABLE

    def _serialize_param_pandas(self, value: Any) -> str | object:
        import pandas as pd

        if isinstance(value, (pd.DataFrame, pd.Series, pd.Index)):
            return repr(pd.util.hash_pandas_object(value).tolist())
        return self._serialize_param_builtin(value)

    @property
    @lru_cache()
    def _serialize_param(self) -> Callable[[Any], str | object]:
        import importlib.util

        if importlib.util.find_spec("pandas") is not None:
            return self._serialize_param_pandas
        else:
            return self._serialize_param_builtin

    @lru_cache()
    def _params(self, job) -> list[Any]:
        return sorted(
            filter(
                lambda p: p is not UNREPRESENTABLE,
                (self._serialize_param(value) for value in job.non_derived_params),
            )
        )

    def all_outputfiles(self):
        from snakemake.jobs import jobfiles

        return jobfiles(self.dag.needrun_jobs(), "output")

    def all_inputfiles(self):
        from snakemake.jobs import jobfiles

        return jobfiles(self.dag.jobs, "input")

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
        try:
            os.remove(self._iocache_filename)
        except FileNotFoundError:
            pass

    @contextmanager
    def noop(self, *args):
        yield
