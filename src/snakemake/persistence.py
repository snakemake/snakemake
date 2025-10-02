__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import asyncio
from dataclasses import dataclass, field
import hashlib
import os
import shutil
import json
import stat
import tempfile
import time
import pickle
from base64 import urlsafe_b64encode, b64encode
from functools import lru_cache
from itertools import count
from pathlib import Path
from contextlib import contextmanager
from typing import Any, Optional, Set, TYPE_CHECKING

from snakemake_interface_executor_plugins.persistence import (
    PersistenceExecutorInterface,
)
from snakemake_interface_executor_plugins.settings import ExecMode

from snakemake.common.tbdstring import TBDString
import snakemake.exceptions
from snakemake.logging import logger
from snakemake.jobs import jobfiles, Job
from snakemake.utils import listfiles
from snakemake.io import _IOFile, is_flagged, get_flag_value, IOCache
from snakemake_interface_common.exceptions import WorkflowError
from snakemake.settings.types import DeploymentMethod

if TYPE_CHECKING:
    from snakemake.dag import DAG

import lmdb

UNREPRESENTABLE = object()
RECORD_FORMAT_VERSION = 6


class Persistence(PersistenceExecutorInterface):
    def __init__(
        self,
        nolock=False,
        dag=None,
        conda_prefix=None,
        singularity_prefix=None,
        shadow_prefix=None,
        warn_only=False,
        path: Path | None = None,
    ):
        import importlib.util

        self._serialize_param = (
            self._serialize_param_pandas
            if importlib.util.find_spec("pandas") is not None
            else self._serialize_param_builtin
        )

        self._max_len = None

        if path is None:
            self._path = Path(os.path.abspath(".snakemake"))
        else:
            self._path = path
        os.makedirs(self.path, exist_ok=True)

        self._lockdir = os.path.join(self.path, "locks")
        os.makedirs(self._lockdir, exist_ok=True)

        self.dag = dag
        self._lockfile = dict()

        self._metadata_path = os.path.join(self.path, "metadata")
        self._incomplete_path = os.path.join(self.path, "incomplete")

        self.conda_env_archive_path = os.path.join(self.path, "conda-archive")
        self.benchmark_path = os.path.join(self.path, "benchmarks")

        self.source_cache = os.path.join(self.path, "source_cache")

        self.iocache_path = os.path.join(self.path, "iocache")

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

        # place to store any auxiliary information needed during a run (e.g. source tarballs)
        self._aux_path = os.path.join(self.path, "auxiliary")

        # migration of .snakemake folder structure
        migration_indicator = Path(
            os.path.join(self._incomplete_path, "migration_underway")
        )
        if (
            os.path.exists(self._metadata_path)
            and not os.path.exists(self._incomplete_path)
        ) or migration_indicator.exists():
            os.makedirs(self._incomplete_path, exist_ok=True)

            migration_indicator.touch()

            self.migrate_v1_to_v2()

            migration_indicator.unlink()

        self._incomplete_cache = None

        for d in (
            self._metadata_path,
            self._incomplete_path,
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

        self._read_record = self._read_record_cached
        self.max_checksum_file_size = (
            self.dag.workflow.dag_settings.max_checksum_file_size
        )

    @property
    def path(self) -> Path:
        return Path(self._path)

    @property
    def aux_path(self) -> Path:
        return Path(self._aux_path)

    def migrate_v1_to_v2(self):
        logger.info("Migrating .snakemake folder to new format...")
        i = 0
        for path, _, filenames in os.walk(self._metadata_path):
            path = Path(path)
            for filename in filenames:
                with open(path / filename, "r") as f:
                    try:
                        record = json.load(f)
                    except json.JSONDecodeError:
                        continue  # not a properly formatted JSON file

                    if record.get("incomplete", False):
                        target_path = Path(self._incomplete_path) / path.relative_to(
                            self._metadata_path
                        )
                        os.makedirs(target_path, exist_ok=True)
                        shutil.copyfile(path / filename, target_path / filename)
                i += 1
                # this can take a while for large folders...
                if (i % 10000) == 0 and i > 0:
                    logger.info(f"{i} files migrated")

        logger.info("Migration complete")

    @property
    def files(self):
        if self._files is None:
            self._files = set(self.dag.output_files)
        return self._files

    @property
    def locked(self):
        inputfiles = set(self.all_inputfiles())
        outputfiles = set(self.all_outputfiles())
        if os.path.exists(self._lockdir):
            for lockfile in self._locks("input"):
                with open(lockfile) as lock:
                    for f in lock:
                        f = f.strip()
                        if f in outputfiles:
                            return True
            for lockfile in self._locks("output"):
                with open(lockfile) as lock:
                    for f in lock:
                        f = f.strip()
                        if f in outputfiles or f in inputfiles:
                            return True
        return False

    @contextmanager
    def lock_warn_only(self):
        if self.locked:
            logger.info(
                "Error: Directory cannot be locked. This usually "
                "means that another Snakemake instance is running on this directory. "
                "Another possibility is that a previous run exited unexpectedly."
            )
        yield

    @contextmanager
    def lock(self):
        if self.locked:
            raise snakemake.exceptions.LockException()
        try:
            self._lock(self.all_inputfiles(), "input")
            self._lock(self.all_outputfiles(), "output")
            yield
        finally:
            self.unlock()

    def unlock(self):
        logger.debug("unlocking")
        for lockfile in self._lockfile.values():
            try:
                logger.debug("removing lock")
                os.remove(lockfile)
            except OSError as e:
                if e.errno != 2:  # missing file
                    raise e
        logger.debug("removed all locks")

    def cleanup_locks(self):
        shutil.rmtree(self._lockdir)

    def cleanup_metadata(self, path):
        return self._delete_record(self._incomplete_path, path) or self._delete_record(
            self._metadata_path, path
        )

    def cleanup_shadow(self):
        if os.path.exists(self.shadow_path):
            shutil.rmtree(self.shadow_path)
            os.mkdir(self.shadow_path)

    def cleanup_containers(self):
        from humanfriendly import format_size

        required_imgs = {Path(img.path) for img in self.dag.container_imgs.values()}
        img_dir = Path(self.container_img_path)
        total_size_cleaned_up = 0
        num_containers_removed = 0
        for pulled_img in img_dir.glob("*.simg"):
            if pulled_img in required_imgs:
                continue
            size_bytes = pulled_img.stat().st_size
            total_size_cleaned_up += size_bytes
            filesize = format_size(size_bytes)
            pulled_img.unlink()
            logger.debug(f"Removed unrequired container {pulled_img} ({filesize})")
            num_containers_removed += 1

        if num_containers_removed == 0:
            logger.info("No containers require cleaning up")
        else:
            logger.info(
                f"Cleaned up {num_containers_removed} containers, saving {format_size(total_size_cleaned_up)}"
            )

    def conda_cleanup_envs(self):
        # cleanup envs
        for address in set(
            env.address
            for env in self.dag.conda_envs.values()
            if not env.is_externally_managed
        ):
            removed = False
            if os.path.exists(address):
                try:
                    shutil.rmtree(address)
                except Exception as e:
                    raise WorkflowError(f"Failed to remove conda env {address}: {e}")
                removed = True
            yaml_path = Path(address).with_suffix(".yaml")
            if yaml_path.exists():
                try:
                    yaml_path.unlink()
                except Exception as e:
                    raise WorkflowError(
                        f"Failed to remove conda env yaml {yaml_path}: {e}"
                    )

                removed = True
            if removed:
                logger.info(f"Removed conda env {address}")

        # cleanup env archives
        in_use = set(env.content_hash for env in self.dag.conda_envs.values())
        for d in os.listdir(self.conda_env_archive_path):
            if d not in in_use:
                shutil.rmtree(os.path.join(self.conda_env_archive_path, d))

    def started(self, job, external_jobid: Optional[str] = None):
        for f in job.output:
            self._record(self._incomplete_path, {"external_jobid": external_jobid}, f)

    def _remove_incomplete_marker(self, job):
        for f in job.output:
            self._delete_record(self._incomplete_path, f)

    async def finished(self, job):
        if not self.dag.workflow.execution_settings.keep_metadata:
            self._remove_incomplete_marker(job)
            # do not store metadata if not requested
            return

        if (
            self.dag.workflow.exec_mode == ExecMode.DEFAULT
            or self.dag.workflow.remote_execution_settings.immediate_submit
        ):
            code = self._code(job.rule)
            input = self._input(job)
            log = self._log(job)
            params = self._params(job)
            shellcmd = job.shellcmd
            conda_env = self._conda_env(job)
            software_stack_hash = self._software_stack_hash(job)
            fallback_time = time.time()
            for f in job.output:
                rec_path = self._record_path(self._incomplete_path, f)
                starttime = (
                    os.path.getmtime(rec_path) if os.path.exists(rec_path) else None
                )
                # Sometimes finished is called twice, if so, lookup the previous starttime
                if not os.path.exists(rec_path):
                    starttime = self._read_record(self._metadata_path, f).get(
                        "starttime", None
                    )

                endtime = (
                    (await f.mtime()).local_or_storage()
                    if await f.exists()
                    else fallback_time
                )

                checksums = (
                    (infile, await infile.checksum(self.max_checksum_file_size))
                    for infile in job.input
                )
                self._record(
                    self._metadata_path,
                    {
                        "record_format_version": RECORD_FORMAT_VERSION,
                        "code": code,
                        "rule": job.rule.name,
                        "input": input,
                        "log": log,
                        "params": params,
                        "shellcmd": shellcmd,
                        "incomplete": False,
                        "starttime": starttime,
                        "endtime": endtime,
                        "job_hash": hash(job),
                        "conda_env": conda_env,
                        "software_stack_hash": software_stack_hash,
                        "container_img_url": job.container_img_url,
                        "input_checksums": {
                            infile: checksum
                            async for infile, checksum in checksums
                            if checksum is not None
                        },
                    },
                    f,
                )
        # remove incomplete marker only after creation of metadata record.
        # otherwise the job starttime will be missing.
        self._remove_incomplete_marker(job)

    def cleanup(self, job):
        for f in job.output:
            self.cleanup_metadata(f)

    async def incomplete(self, job):
        if self._incomplete_cache is None:
            self._cache_incomplete_folder()

        if self._incomplete_cache is False:  # cache deactivated

            def marked_incomplete(f):
                return self._exists_record(self._incomplete_path, f)

        else:

            def marked_incomplete(f):
                rec_path = self._record_path(self._incomplete_path, f)
                return rec_path in self._incomplete_cache

        async def is_incomplete(f):
            exists = await f.exists()
            marked = marked_incomplete(f)
            return f if exists and marked else None

        async with asyncio.TaskGroup() as tg:
            tasks = [tg.create_task(is_incomplete(f)) for f in job.output]

        return [task.result() for task in tasks]

    def _cache_incomplete_folder(self):
        self._incomplete_cache = {
            os.path.join(path, f)
            for path, dirnames, filenames in os.walk(self._incomplete_path)
            for f in filenames
        }

    def external_jobids(self, job):
        return list(
            set(
                self._read_record(self._incomplete_path, f).get("external_jobid", None)
                for f in job.output
            )
        )

    def has_metadata(self, job: Job) -> bool:
        return all(self.metadata(path) for path in job.output)

    def has_outdated_metadata(self, job: Job) -> bool:
        return any(
            self.metadata(path).get("record_format_version", 0) < RECORD_FORMAT_VERSION
            for path in job.output
        )

    def metadata(self, path):
        return self._read_record(self._metadata_path, path)

    def rule(self, path):
        return self.metadata(path).get("rule")

    def input(self, path):
        return self.metadata(path).get("input")

    def log(self, path):
        return self.metadata(path).get("log")

    def shellcmd(self, path):
        return self.metadata(path).get("shellcmd")

    def params(self, path):
        return self.metadata(path).get("params")

    def code(self, path):
        return self.metadata(path).get("code")

    def record_format_version(self, path):
        return self.metadata(path).get("record_format_version")

    def conda_env(self, path):
        return self.metadata(path).get("conda_env")

    def container_img_url(self, path):
        return self.metadata(path).get("container_img_url")

    def input_checksums(self, job, input_path):
        """Return all checksums of the given input file
        recorded for the output of the given job.
        """
        return set(
            self.metadata(output_path).get("input_checksums", {}).get(input_path)
            for output_path in job.output
        )

    def code_changed(self, job, file=None):
        """Yields output files with changed code or bool if file given."""
        return _bool_or_gen(self._code_changed, job, file=file)

    def input_changed(self, job, file=None):
        """Yields output files with changed input or bool if file given."""
        return _bool_or_gen(self._input_changed, job, file=file)

    def params_changed(self, job, file=None):
        """Yields output files with changed params or bool if file given."""
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

    def software_stack_changed(self, job, file=None):
        """Yields output files with changed software env or bool if file given."""
        return _bool_or_gen(self._software_stack_changed, job, file=file)

    def _code_changed(self, job, file=None):
        assert file is not None
        fmt_version = self.record_format_version(file)
        if fmt_version is None or fmt_version < 3:
            # no reliable code stored
            return False
        recorded = self.code(file)
        return recorded is not None and recorded != self._code(job.rule)

    def _input_changed(self, job, file=None):
        assert file is not None
        fmt_version = self.record_format_version(file)
        if fmt_version is None or fmt_version < 4:
            # no reliable input stored
            return False
        recorded = self.input(file)
        return recorded is not None and recorded != self._input(job)

    def _software_stack_changed(self, job, file=None):
        assert file is not None
        fmt_version = self.record_format_version(file)
        if fmt_version is None or fmt_version < 5:
            # no reliable software stack hash stored (previous storage ignored pin files
            # and aux deploy files of conda envs as well as env modules)
            return False

        recorded = self.software_stack_hash(file)
        return recorded is not None and recorded != self._software_stack_hash(job)

    def software_stack_hash(self, path):
        return self.metadata(path).get("software_stack_hash")

    def _software_stack_hash(self, job):
        # TODO move code for retrieval into software deployment plugin interface once
        # available
        md5hash = hashlib.md5(usedforsecurity=False)
        if (
            DeploymentMethod.CONDA
            in self.dag.workflow.deployment_settings.deployment_method
            and job.conda_env
        ):
            md5hash.update(job.conda_env.hash.encode())
        if (
            DeploymentMethod.APPTAINER
            in self.dag.workflow.deployment_settings.deployment_method
            and job.container_img_url
        ):
            md5hash.update(job.container_img_url.encode())
        if job.env_modules:
            md5hash.update(job.env_modules.hash.encode())
        return md5hash.hexdigest()

    @contextmanager
    def noop(self, *args):
        yield

    def _b64id(self, s):
        return urlsafe_b64encode(str(s).encode()).decode()

    @lru_cache()
    def _code(self, rule):
        # Scripts and notebooks are triggered by changes in the script mtime.
        # Changes to python and shell rules are triggered by changes in the plain text.
        if rule.shellcmd is not None:
            return rule.shellcmd
        if rule.run_func_src is not None:
            return rule.run_func_src
        return None

    @lru_cache()
    def _conda_env(self, job):
        if job.conda_env:
            return b64encode(job.conda_env.content).decode()

    @lru_cache()
    def _input(self, job):
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
    def _log(self, job):
        return sorted(job.log)

    def _serialize_param_builtin(self, value: Any):
        if (
            value is None
            or isinstance(
                value,
                (
                    int,
                    float,
                    bool,
                    str,
                    complex,
                    range,
                    list,
                    tuple,
                    dict,
                    set,
                    frozenset,
                    bytes,
                    bytearray,
                ),
            )
            and value is not TBDString
        ):
            return repr(value)
        else:
            return UNREPRESENTABLE

    def _serialize_param_pandas(self, value: Any):
        import pandas as pd

        if isinstance(value, (pd.DataFrame, pd.Series, pd.Index)):
            return repr(pd.util.hash_pandas_object(value).tolist())
        return self._serialize_param_builtin(value)

    @lru_cache()
    def _params(self, job: Job):
        return sorted(
            filter(
                lambda p: p is not UNREPRESENTABLE,
                (self._serialize_param(value) for value in job.non_derived_params),
            )
        )

    @lru_cache()
    def _output(self, job):
        return sorted(job.output)

    def _record(
        self,
        subject,
        json_value,
        id,
        mode=stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP,
    ):
        recpath = self._record_path(subject, id)
        try:
            recpath_stat = os.stat(recpath)
        except FileNotFoundError:
            recpath_stat = None
            recdir = os.path.dirname(recpath)
            os.makedirs(recdir, exist_ok=True)

        with open(recpath, "w") as recfile:
            json.dump(json_value, recfile)

        # ensure read and write permissions for user and group if they don't include the required mode
        if recpath_stat is None:
            os.chmod(recpath, mode)
        else:
            existing = stat.S_IMODE(recpath_stat.st_mode)
            new_mode = existing | mode
            if existing != new_mode:
                os.chmod(recpath, new_mode)

    def _delete_record(self, subject, id):
        try:
            recpath = self._record_path(subject, id)
            os.remove(recpath)
            recdirs = os.path.relpath(os.path.dirname(recpath), start=subject)
            if recdirs != ".":
                os.removedirs(recdirs)
            return True
        except OSError as e:
            if e.errno != 2:
                # not missing
                raise e
            else:
                # file is missing, report failure
                return False

    @lru_cache()
    def _read_record_cached(self, subject, id):
        return self._read_record_uncached(subject, id)

    def _read_record_uncached(self, subject, id):
        if not self._exists_record(subject, id):
            return dict()
        path = self._record_path(subject, id)
        with open(path, "r") as f:
            try:
                return json.load(f)
            except json.JSONDecodeError:
                # Since record writing cannot be reliably made atomic (some network
                # filesystems, e.g. gluster have issues with writing to a temp file
                # and then moving) we ignore corrupted or incompletely written records
                # here.
                # They can only occur if a snakemake process is running and one does a
                # dry-run (or intentionally disables locking) at the same time.
                logger.warning(
                    f"Ignore corrupted or currently written metadata record {path}."
                )
                return dict()

    def _exists_record(self, subject, id):
        return os.path.exists(self._record_path(subject, id))

    def _locks(self, type):
        return (
            f
            for f, _ in listfiles(
                os.path.join(self._lockdir, f"{{n,[0-9]+}}.{type}.lock")
            )
            if not os.path.isdir(f)
        )

    def _lock(self, files, type):
        for i in count(0):
            lockfile = os.path.join(self._lockdir, f"{i}.{type}.lock")
            if not os.path.exists(lockfile):
                self._lockfile[type] = lockfile
                with open(lockfile, "w") as lock:
                    print(*files, sep="\n", file=lock)
                return

    def _fetch_max_len(self, subject):
        if self._max_len is None:
            self._max_len = os.pathconf(subject, "PC_NAME_MAX")
        return self._max_len

    def _record_path(self, subject, id: _IOFile):
        assert isinstance(id, _IOFile)
        id = id.storage_object.query if id.is_storage else id

        max_len = (
            self._fetch_max_len(subject) if os.name == "posix" else 255
        )  # maximum NTFS and FAT32 filename length
        if max_len == 0:
            max_len = 255

        b64id = self._b64id(id)
        # split into chunks of proper length
        b64id = [b64id[i : i + max_len - 1] for i in range(0, len(b64id), max_len - 1)]
        # prepend dirs with @ (does not occur in b64) to avoid conflict with b64-named files in the same dir
        b64id = ["@" + s for s in b64id[:-1]] + [b64id[-1]]
        path = os.path.join(subject, *b64id)
        return path

    def all_outputfiles(self):
        # we only look at output files that will be updated
        return jobfiles(self.dag.needrun_jobs(), "output")

    def all_inputfiles(self):
        # we consider all input files, also of not running jobs
        return jobfiles(self.dag.jobs, "input")

    def deactivate_cache(self):
        self._read_record_cached.cache_clear()
        self._read_record = self._read_record_uncached
        self._incomplete_cache = False

    @property
    def _iocache_filename(self):
        return os.path.join(self.iocache_path, "latest.pkl")

    def save_iocache(self):
        filepath = self._iocache_filename
        with open(filepath, "wb") as handle:
            self.dag.workflow.iocache.save(handle)

    def load_iocache(self):
        filepath = self._iocache_filename
        if os.path.exists(filepath):
            logger.info("Loading trusted IOCache from latest dry-run.")
            with open(filepath, "rb") as handle:
                self.dag.workflow.iocache = IOCache.load(handle)

    def drop_iocache(self):
        filepath = self._iocache_filename
        if os.path.exists(filepath):
            os.remove(filepath)


def _bool_or_gen(func, job, file=None):
    if file is None:
        return (f for f in job.output if func(job, file=f))
    else:
        return func(job, file=file)


@dataclass
class ParamsChange:
    only_old: Set[Any] = field(default_factory=set)
    only_new: Set[Any] = field(default_factory=set)
    files: Set[str] = field(default_factory=set)

    def __post_init__(self):
        if not self:
            self.files = set()

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


class LmdbPersistence(Persistence):
    """
    LMDB-based persistence implementation for Snakemake metadata.

    Uses Lightning Memory-Mapped Database (LMDB) for efficient storage and retrieval
    of job metadata, incomplete markers, and other persistence data.
    """

    def __init__(
        self,
        dag: "DAG",
        nolock=False,
        conda_prefix=None,
        singularity_prefix=None,
        shadow_prefix=None,
        warn_only=False,
        path: Path | None = None,
    ):
        import importlib.util

        self._serialize_param = (
            self._serialize_param_pandas
            if importlib.util.find_spec("pandas") is not None
            else self._serialize_param_builtin
        )

        self.dag = dag

        if path is None:
            self._path = Path(os.path.abspath(".snakemake"))
        else:
            self._path = path
        os.makedirs(self.path, exist_ok=True)

        self._lmdb_path = os.path.join(self.path, "lmdb")

        self.conda_env_archive_path = os.path.join(self.path, "conda-archive")
        self.benchmark_path = os.path.join(self.path, "benchmarks")
        self.source_cache = os.path.join(self.path, "source_cache")
        self.iocache_path = os.path.join(self.path, "iocache")

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

        self._aux_path = os.path.join(self.path, "auxiliary")

        for d in (
            self.shadow_path,
            self.conda_env_archive_path,
            self.conda_env_path,
            self.container_img_path,
            self.aux_path,
            self.iocache_path,
            self._lmdb_path,
        ):
            os.makedirs(d, exist_ok=True)

        self._env = self._init_lmdb()

        # Lock handling
        self._locked = False
        if nolock:
            self.lock = self.noop
            self.unlock = self.noop_unlock
        if warn_only:
            self.lock = self.lock_warn_only
            self.unlock = self.noop_unlock

        self._incomplete_cache = None
        self.max_checksum_file_size = (
            self.dag.workflow.dag_settings.max_checksum_file_size
        )

    @property
    def path(self) -> Path:
        return Path(self._path)

    @property
    def aux_path(self) -> Path:
        return Path(self._aux_path)

    @contextmanager
    def noop(self, *args):
        """No-op context manager for disabled locking."""
        yield

    def noop_unlock(self) -> None:
        """No-op unlock for disabled locking."""
        pass

    @contextmanager
    def lock_warn_only(self):
        """Warning-only lock context manager."""
        if self.locked:
            logger.info(
                "Error: Directory cannot be locked. This usually "
                "means that another Snakemake instance is running on this directory. "
                "Another possibility is that a previous run exited unexpectedly."
            )
        yield

    def _init_lmdb(self) -> lmdb.Environment:
        """Initialize LMDB environment and databases."""
        _env: lmdb.Environment = lmdb.open(
            self._lmdb_path,
            max_dbs=4,
            map_size=1024**3,
        )

        with _env.begin(write=True) as txn:
            self._metadata_db = _env.open_db(b"metadata", txn=txn)
            self._incomplete_db = _env.open_db(b"incomplete", txn=txn)
            self._locks_db = _env.open_db(b"locks", txn=txn)
            self._cache_db = _env.open_db(b"cache", txn=txn)
        return _env

    @property
    def locked(self) -> bool:
        if self._locked:
            return True

        inputfiles = set(str(f) for f in self.all_inputfiles())
        outputfiles = set(str(f) for f in self.all_outputfiles())

        try:
            with self._env.begin(db=self._locks_db) as txn:
                cursor = txn.cursor()
                for key, value in cursor:
                    if key in [b"input_lock", b"output_lock"]:
                        locked_files = set(json.loads(value.decode("utf-8")))
                        if key == b"output_lock":
                            if locked_files & (outputfiles | inputfiles):
                                return True
                        elif key == b"input_lock":
                            if locked_files & outputfiles:
                                return True
        except Exception:
            # If LMDB access fails, fall back to parent implementation
            return super().locked
        return False

    @contextmanager
    def lock(self):
        """Context manager for locking."""
        if self.locked:
            raise snakemake.exceptions.LockException()
        try:
            self._lock_files(self.all_inputfiles(), "input")
            self._lock_files(self.all_outputfiles(), "output")
            yield
        finally:
            self.unlock()

    def _lock_files(self, files, lock_type: str):
        """Store lock information in LMDB."""
        with self._env.begin(write=True) as txn:
            lock_key = f"{lock_type}_lock".encode("utf-8")
            files_list = list(str(f) for f in files)
            lock_data = json.dumps(files_list).encode("utf-8")
            txn.put(lock_key, lock_data, db=self._locks_db)
        self._locked = True

    def unlock(self) -> None:
        """Release locks."""
        self._locked = False
        with self._env.begin(write=True) as txn:
            txn.drop(db=self._locks_db, delete=False)

    @lru_cache()
    def _file_key(self, file: _IOFile) -> bytes:
        """Generate LMDB key for a file."""
        file_str = file.storage_object.query if file.is_storage else str(file)
        return self._b64id(file_str).encode("utf-8")

    def _store_record(self, db, record_data: dict[str, Any], file: _IOFile) -> None:
        """Store a record in the specified LMDB database."""
        if not (
            self.dag.workflow.exec_mode == ExecMode.DEFAULT
            or self.dag.workflow.remote_execution_settings.immediate_submit
        ):
            return
        key = self._file_key(file)
        value = pickle.dumps(record_data)
        with self._env.begin(write=True) as txn:
            txn.put(key, value, db=db)

    def _get_record(self, db, file: _IOFile) -> dict[str, Any]:
        """Get a record from the specified LMDB database."""
        key = self._file_key(file)
        with self._env.begin(db=db) as txn:
            value = txn.get(key)
            if value is None:
                return {}
            try:
                return pickle.loads(value)
            except (pickle.UnpicklingError, AttributeError):
                logger.warning(f"Corrupted metadata record for {file}")
                return {}

    def cleanup_locks(self) -> None:
        """Clean up all locks."""
        with self._env.begin(write=True) as txn:
            txn.drop(db=self._locks_db, delete=False)

    def cleanup_metadata(self, path: _IOFile) -> bool:
        """Remove metadata for a file."""
        key = self._file_key(path)
        with self._env.begin(write=True) as txn:
            deleted_metadata = txn.delete(key, db=self._metadata_db)
            deleted_incomplete = txn.delete(key, db=self._incomplete_db)
            return deleted_metadata or deleted_incomplete

    def started(self, job: Job, external_jobid: str | None = None) -> None:
        """Mark job as started by storing incomplete markers."""

        for f in job.output:
            self._store_record(
                self._incomplete_db, {"external_jobid": external_jobid}, f
            )

    async def finished(self, job: Job) -> None:
        """Mark job as finished and store metadata."""
        if not self.dag.workflow.execution_settings.keep_metadata:
            self._remove_incomplete_marker(job)
            return

        if (
            self.dag.workflow.exec_mode == ExecMode.DEFAULT
            or self.dag.workflow.remote_execution_settings.immediate_submit
        ):
            code = self._code(job.rule)
            input = self._input(job)
            log = self._log(job)
            params = self._params(job)
            shellcmd = job.shellcmd
            conda_env = self._conda_env(job)
            software_stack_hash = self._software_stack_hash(job)
            fallback_time = time.time()

            # here we just make a new transcation so that all records are in 1 txn for this job
            # rather than using _store_recordm

            with self._env.begin(write=True) as txn:
                for f in job.output:
                    starttime = self._get_starttime(f, fallback_time)
                    endtime = (
                        (await f.mtime()).local_or_storage()
                        if await f.exists()
                        else fallback_time
                    )

                    checksums = (
                        (infile, await infile.checksum(self.max_checksum_file_size))
                        for infile in job.input
                    )

                    record_data = {
                        "record_format_version": RECORD_FORMAT_VERSION,
                        "code": code,
                        "rule": job.rule.name,
                        "input": input,
                        "log": log,
                        "params": params,
                        "shellcmd": shellcmd,
                        "incomplete": False,
                        "starttime": starttime,
                        "endtime": endtime,
                        "job_hash": hash(job),
                        "conda_env": conda_env,
                        "software_stack_hash": software_stack_hash,
                        "container_img_url": job.container_img_url,
                        "input_checksums": {
                            str(infile): checksum
                            async for infile, checksum in checksums
                            if checksum is not None
                        },
                    }

                    key = self._file_key(f)
                    value = pickle.dumps(record_data)
                    txn.put(key, value, db=self._metadata_db)

        self._remove_incomplete_marker(job)

    def _get_starttime(self, file: _IOFile, fallback_time: float) -> float | None:
        """Get start time from incomplete marker or existing metadata."""
        # First try to get from incomplete marker
        incomplete_record = self._get_record(self._incomplete_db, file)
        if incomplete_record and "starttime" in incomplete_record:
            return incomplete_record["starttime"]

        # If not found, try existing metadata
        metadata_record = self._get_record(self._metadata_db, file)
        if metadata_record and "starttime" in metadata_record:
            return metadata_record["starttime"]

        return fallback_time

    def _remove_incomplete_marker(self, job: Job) -> None:
        """Remove incomplete markers for job outputs."""
        for f in job.output:
            key = self._file_key(f)
            with self._env.begin(write=True) as txn:
                txn.delete(key, db=self._incomplete_db)

    async def incomplete(self, job: Job) -> list[_IOFile | None]:
        """Return list of incomplete output files for the job."""
        if self._incomplete_cache is None:
            self._cache_incomplete_files()

        async def is_incomplete(f):
            exists = await f.exists()
            marked = self._is_marked_incomplete(f)
            return f if exists and marked else None

        async with asyncio.TaskGroup() as tg:
            tasks = [tg.create_task(is_incomplete(f)) for f in job.output]

        return [task.result() for task in tasks]

    def _is_marked_incomplete(self, file: _IOFile) -> bool:
        """Check if file is marked as incomplete."""
        if self._incomplete_cache is False:  # cache deactivated
            key = self._file_key(file)
            with self._env.begin(db=self._incomplete_db) as txn:
                return txn.get(key) is not None
        else:
            key = self._file_key(file)
            return key in self._incomplete_cache

    def _cache_incomplete_files(self) -> None:
        """Cache all incomplete file keys for performance."""
        self._incomplete_cache = set()
        with self._env.begin(db=self._incomplete_db) as txn:
            cursor = txn.cursor()
            for key, _ in cursor:
                self._incomplete_cache.add(key)

    @lru_cache()
    def _input(self, job):
        """Override to ensure all inputs are converted to strings for pickling."""

        def get_paths():
            for f in job.input:
                if f.is_storage:
                    yield f.storage_object.query
                elif is_flagged(f, "pipe"):
                    yield "<pipe>"
                elif is_flagged(f, "service"):
                    yield "<service>"
                else:
                    # Convert to string to avoid pickling _IOFile objects with callable attributes
                    path = (
                        get_flag_value(f, "sourcecache_entry")
                        if is_flagged(f, "sourcecache_entry")
                        else f
                    )
                    yield str(path)

        return sorted(get_paths())

    @lru_cache()
    def _log(self, job):
        """Override to ensure all log files are converted to strings for pickling."""
        return sorted(str(f) for f in job.log)

    # Note: external_jobids is inherited from parent class - it can use our overridden metadata() method

    def metadata(self, path: _IOFile) -> dict[str, Any]:
        return self._get_record(self._metadata_db, path)

    # Note: input_checksums, _code_changed, _input_changed, _software_stack_changed,
    # and deactivate_cache are inherited from parent class since they use the same logic

    def __del__(self):
        """Clean up LMDB environment on deletion."""
        if hasattr(self, "_env") and self._env:
            self._env.close()
