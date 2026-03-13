__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import shutil
import json
import stat
from base64 import urlsafe_b64encode
from functools import lru_cache
from itertools import count
from pathlib import Path
from typing import Iterable

from snakemake.persistence import PersistenceBase, MetadataRecord
from snakemake.utils import listfiles
from snakemake.logging import logger


class FilePersistence(PersistenceBase):
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
        super().__init__(
            nolock=nolock,
            dag=dag,
            conda_prefix=conda_prefix,
            singularity_prefix=singularity_prefix,
            shadow_prefix=shadow_prefix,
            warn_only=warn_only,
            path=path,
        )
        self._max_len = None
        self._incomplete_cache = None

        self._lockdir = os.path.join(self.path, "locks")
        os.makedirs(self._lockdir, exist_ok=True)
        self._lockfile = dict()

        self._metadata_path = os.path.join(self.path, "metadata")
        self._incomplete_path = os.path.join(self.path, "incomplete")

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

        for d in (self._metadata_path, self._incomplete_path):
            os.makedirs(d, exist_ok=True)

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

    def _fetch_max_len(self, subject):
        if self._max_len is None:
            self._max_len = os.pathconf(subject, "PC_NAME_MAX")
        return self._max_len

    def _record_path(self, subject, id: str):
        max_len = self._fetch_max_len(subject) if os.name == "posix" else 255
        if max_len == 0:
            max_len = 255
        b64id = urlsafe_b64encode(id.encode()).decode()
        b64id = [b64id[i : i + max_len - 1] for i in range(0, len(b64id), max_len - 1)]
        b64id = ["@" + s for s in b64id[:-1]] + [b64id[-1]]
        return os.path.join(subject, *b64id)

    def _io_write(
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

    def _io_delete(self, subject, id):
        try:
            recpath = self._record_path(subject, id)
            os.remove(recpath)

            # also remove the directory if it is empty to avoid leaving a large number of empty directories around
            recdir = os.path.dirname(recpath)
            if recdir != subject:
                try:
                    os.removedirs(recdir)
                except OSError:
                    # ignore errors (e.g., when the directory is not empty or already removed by another process)
                    pass
            return True
        except OSError as e:
            if e.errno != 2:
                # not missing
                raise e
            else:
                # file is missing, report failure
                return False

    def _io_read(self, subject, id):
        path = self._record_path(subject, id)
        if not os.path.exists(path):
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

    def _clear_cache(self) -> None:
        self._read_record_cached.cache_clear()
        self._incomplete_cache = None

    @lru_cache()
    def _read_record_cached(self, key: str) -> MetadataRecord | None:
        rec = self._io_read(self._metadata_path, key)
        return MetadataRecord.model_validate(rec) if rec else None

    def _read_record(self, key: str) -> MetadataRecord | None:
        return self._read_record_cached(key)

    def _write_record(self, key: str, record: MetadataRecord) -> None:
        self._io_write(self._metadata_path, dict(record.items()), key)

    def _delete_record(self, key: str) -> bool:
        return self._io_delete(self._metadata_path, key)

    def _mark_incomplete(self, key: str, external_jobid: str | None) -> None:
        self._io_write(self._incomplete_path, {"external_jobid": external_jobid}, key)
        if self._incomplete_cache is not None:
            self._incomplete_cache.add(self._record_path(self._incomplete_path, key))

    def _unmark_incomplete(self, key: str) -> None:
        self._io_delete(self._incomplete_path, key)
        if self._incomplete_cache is not None:
            self._incomplete_cache.discard(
                self._record_path(self._incomplete_path, key)
            )

    def _filter_incomplete_keys(self, keys: Iterable[str]) -> set[str]:
        if self._incomplete_cache is None:
            self._incomplete_cache = {
                os.path.join(path, f)
                for path, _, filenames in os.walk(self._incomplete_path)
                for f in filenames
            }

        return {
            k
            for k in keys
            if self._record_path(self._incomplete_path, k) in self._incomplete_cache
        }

    def _get_external_jobids(self, keys: Iterable[str]) -> set[str]:
        jobids = set()
        for k in keys:
            rec = self._io_read(self._incomplete_path, k)
            if "external_jobid" in rec and rec["external_jobid"] is not None:
                jobids.add(rec["external_jobid"])
        return jobids

    def _read_locks(self) -> Iterable[tuple[str, str]]:
        for lock_type in ["input", "output"]:
            for lockfile, _ in listfiles(
                os.path.join(self._lockdir, f"{{n,[0-9]+}}.{lock_type}.lock")
            ):
                if not os.path.isdir(lockfile):
                    with open(lockfile) as lock:
                        for f in lock:
                            yield lock_type, f.strip()

    def _write_locks(self, lock_type: str, keys: Iterable[str]) -> None:
        keys_list = list(keys)
        if not keys_list:
            return
        for i in count(0):
            lockfile = os.path.join(self._lockdir, f"{i}.{lock_type}.lock")
            if not os.path.exists(lockfile):
                self._lockfile[lock_type] = lockfile
                with open(lockfile, "w") as lock:
                    print(*keys_list, sep="\n", file=lock)
                return

    def _delete_locks(self) -> None:
        for lockfile in self._lockfile.values():
            try:
                os.remove(lockfile)
            except OSError as e:
                if e.errno != 2:
                    raise e
        self._lockfile.clear()
        if os.path.exists(self._lockdir):
            shutil.rmtree(self._lockdir)
