import hashlib
import os
import shutil
from abc import abstractmethod
from base64 import b64encode
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, ContextManager, Dict, List, Optional, Set, Union, Callable

from snakemake_interface_executor_plugins.persistence import (
    PersistenceExecutorInterface,
)

from snakemake.common.tbdstring import TBDString
from snakemake.io import get_flag_value, is_flagged
from snakemake.settings.types import DeploymentMethod
from snakemake_interface_common.exceptions import WorkflowError
from snakemake.logging import logger

RECORD_FORMAT_VERSION = 6
UNREPRESENTABLE = object()


@dataclass
class MetadataRecord:
    rule: Optional[str] = None
    input: Optional[List[str]] = None
    log: Optional[List[str]] = None
    shellcmd: Optional[str] = None
    params: Optional[List[Any]] = None
    code: Optional[str] = None
    record_format_version: int = 0
    conda_env: Optional[str] = None
    container_img_url: Optional[str] = None
    software_stack_hash: Optional[str] = None
    job_hash: Optional[int] = None
    starttime: Optional[float] = None
    endtime: Optional[float] = None
    incomplete: Optional[bool] = None
    external_jobid: Optional[str] = None
    input_checksums: Optional[Dict[str, Any]] = field(default_factory=dict)

    def __getitem__(self, key: str) -> Any:
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(key)

    def get(self, key: str, default: Any = None) -> Any:
        return getattr(self, key, default)

    def keys(self):
        return self.__dict__.keys()

    def items(self):
        return self.__dict__.items()


@dataclass
class ParamsChange:
    only_old: Set[Any] = field(default_factory=set)
    only_new: Set[Any] = field(default_factory=set)
    files: Set[str] = field(default_factory=set)

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

    shadow_path: Union[str, Path]
    container_img_path: Union[str, Path]
    conda_env_archive_path: Union[str, Path]
    dag: Any

    def cleanup_shadow(self) -> None:
        if os.path.exists(self.shadow_path):
            shutil.rmtree(self.shadow_path)
            os.mkdir(self.shadow_path)

    def cleanup_containers(self) -> None:
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

    def conda_cleanup_envs(self) -> None:
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


class AbstractPersistence(
    PersistenceExecutorInterface, EnvironmentMaintenanceMixin, FileBackupMixin
):
    dag: Any

    @abstractmethod
    def metadata(self, path: Any) -> Optional[MetadataRecord]: ...

    @property
    @abstractmethod
    def locked(self) -> bool: ...

    @abstractmethod
    def lock_warn_only(self) -> ContextManager[None]: ...

    @abstractmethod
    def lock(self) -> ContextManager[None]: ...

    @abstractmethod
    def unlock(self) -> None: ...

    @abstractmethod
    def cleanup_locks(self) -> None: ...

    @abstractmethod
    def cleanup_metadata(self, path: Any) -> bool: ...

    @abstractmethod
    def started(self, job: Any, external_jobid: Optional[str] = None) -> None: ...

    @abstractmethod
    async def finished(self, job: Any) -> None: ...

    @abstractmethod
    async def incomplete(self, job: Any) -> List[Any]: ...

    @abstractmethod
    def external_jobids(self, job: Any) -> List[str]: ...

    @abstractmethod
    def deactivate_cache(self) -> None: ...

    @abstractmethod
    def save_iocache(self) -> None: ...

    @abstractmethod
    def load_iocache(self) -> None: ...

    @abstractmethod
    def drop_iocache(self) -> None: ...

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

    def rule(self, path: Any) -> Optional[str]:
        return (m := self.metadata(path)) and m.rule

    def input(self, path: Any) -> Optional[List[str]]:
        return (m := self.metadata(path)) and m.input

    def log(self, path: Any) -> Optional[List[str]]:
        return (m := self.metadata(path)) and m.log

    def shellcmd(self, path: Any) -> Optional[str]:
        return (m := self.metadata(path)) and m.shellcmd

    def params(self, path: Any) -> Optional[List[Any]]:
        return (m := self.metadata(path)) and m.params

    def code(self, path: Any) -> Optional[str]:
        return (m := self.metadata(path)) and m.code

    def record_format_version(self, path: Any) -> Optional[int]:
        return (m := self.metadata(path)) and m.record_format_version

    def conda_env(self, path: Any) -> Optional[str]:
        return (m := self.metadata(path)) and m.conda_env

    def container_img_url(self, path: Any) -> Optional[str]:
        return (m := self.metadata(path)) and m.container_img_url

    def software_stack_hash(self, path: Any) -> Optional[str]:
        return (m := self.metadata(path)) and m.software_stack_hash

    def input_checksums(self, job: Any, input_path: Any) -> Set[Any]:
        return set(
            (
                m.input_checksums.get(input_path)
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
            return False
        recorded = self.software_stack_hash(file)
        return recorded is not None and recorded != self._software_stack_hash(job)

    def params_changed(self, job: Any, file=None) -> Union[bool, ParamsChange]:
        files = [file] if file is not None else job.output
        changes = NO_PARAMS_CHANGE
        new = set(self._params(job))

        for outfile in files:
            fmt_version = self.record_format_version(outfile)
            if fmt_version is None or fmt_version < 6:
                continue
            recorded = self.params(outfile)
            if recorded is not None:
                old = set(recorded)
                changes |= ParamsChange(
                    only_old=old - new, only_new=new - old, files={outfile}
                )
        return changes

    @lru_cache()
    def _code(self, rule) -> Optional[str]:
        if rule.shellcmd is not None:
            return rule.shellcmd
        if rule.run_func_src is not None:
            return rule.run_func_src
        return None

    @lru_cache()
    def _conda_env(self, job) -> Optional[str]:
        if job.conda_env:
            return b64encode(job.conda_env.content).decode()

    @lru_cache()
    def _input(self, job) -> List[str]:
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
                        get_flag_value(f, "sourcecache_entry")
                        if is_flagged(f, "sourcecache_entry")
                        else f
                    )

        return sorted(get_paths())

    @lru_cache()
    def _log(self, job) -> List[str]:
        return sorted(job.log)

    def _serialize_param_builtin(self, value: Any) -> Union[str, object]:
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

    def _serialize_param_pandas(self, value: Any) -> Union[str, object]:
        import pandas as pd

        if isinstance(value, (pd.DataFrame, pd.Series, pd.Index)):
            return repr(pd.util.hash_pandas_object(value).tolist())
        return self._serialize_param_builtin(value)

    @property
    @lru_cache()
    def _serialize_param(self) -> Callable[[Any], Union[str, object]]:
        import importlib.util

        if importlib.util.find_spec("pandas") is not None:
            return self._serialize_param_pandas
        else:
            return self._serialize_param_builtin

    @lru_cache()
    def _params(self, job) -> List[Any]:
        return sorted(
            filter(
                lambda p: p is not UNREPRESENTABLE,
                (self._serialize_param(value) for value in job.non_derived_params),
            )
        )

    def _software_stack_hash(self, job) -> str:
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

    def all_outputfiles(self):
        from snakemake.jobs import jobfiles

        return jobfiles(self.dag.needrun_jobs(), "output")

    def all_inputfiles(self):
        from snakemake.jobs import jobfiles

        return jobfiles(self.dag.jobs, "input")
