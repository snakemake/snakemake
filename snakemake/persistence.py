__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import shutil
import signal
import marshal
import pickle
import json
from base64 import urlsafe_b64encode, b64encode
from functools import lru_cache, partial
from itertools import filterfalse, count

from snakemake.logging import logger
from snakemake.jobs import jobfiles
from snakemake.utils import listfiles


class Persistence:
    def __init__(self,
                 nolock=False,
                 dag=None,
                 conda_prefix=None,
                 singularity_prefix=None,
                 shadow_prefix=None,
                 warn_only=False):
        self.path = os.path.abspath(".snakemake")
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        self._lockdir = os.path.join(self.path, "locks")
        if not os.path.exists(self._lockdir):
            os.mkdir(self._lockdir)

        self.dag = dag
        self._lockfile = dict()

        self._metadata_path = os.path.join(self.path, "metadata")
        self.conda_env_archive_path = os.path.join(self.path, "conda-archive")
        self.benchmark_path = os.path.join(self.path, "benchmarks")

        if conda_prefix is None:
            self.conda_env_path = os.path.join(self.path, "conda")
        else:
            self.conda_env_path = os.path.abspath(
                os.path.expanduser(conda_prefix))
        if singularity_prefix is None:
            self.singularity_img_path = os.path.join(self.path, "singularity")
        else:
            self.singularity_img_path = os.path.abspath(
                os.path.expanduser(singularity_prefix))
        if shadow_prefix is None:
            self.shadow_path = os.path.join(self.path, "shadow")
        else:
            self.shadow_path = os.path.join(shadow_prefix, "shadow")

        for d in (self._metadata_path, self.shadow_path, self.conda_env_archive_path, self.conda_env_path, self.singularity_img_path):
            os.makedirs(d, exist_ok=True)

        if nolock:
            self.lock = self.noop
            self.unlock = self.noop
        if warn_only:
            self.lock = self.lock_warn_only
            self.unlock = self.noop

        self._read_record = self._read_record_cached

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

    def lock_warn_only(self):
        if self.locked:
            logger.info(
                "Error: Directory cannot be locked. This usually "
                "means that another Snakemake instance is running on this directory. "
                "Another possibility is that a previous run exited unexpectedly.")

    def lock(self):
        if self.locked:
            raise IOError("Another snakemake process "
                          "has locked this directory.")
        self._lock(self.all_inputfiles(), "input")
        self._lock(self.all_outputfiles(), "output")

    def unlock(self, *args):
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
        self._delete_record(self._metadata_path, path)

    def cleanup_shadow(self):
        if os.path.exists(self.shadow_path):
            shutil.rmtree(self.shadow_path)
            os.mkdir(self.shadow_path)

    def cleanup_conda(self):
        # cleanup envs
        in_use = set(env.hash[:8] for env in self.dag.conda_envs.values())
        for d in os.listdir(self.conda_env_path):
            if len(d)>=8 and d[:8] not in in_use:
                if os.path.isdir(os.path.join(self.conda_env_path, d)):
                    shutil.rmtree(os.path.join(self.conda_env_path, d))
                else:
                    os.remove(os.path.join(self.conda_env_path, d))

        # cleanup env archives
        in_use = set(env.content_hash for env in self.dag.conda_envs.values())
        for d in os.listdir(self.conda_env_archive_path):
            if d not in in_use:
                shutil.rmtree(os.path.join(self.conda_env_archive_path, d))

    def started(self, job, external_jobid=None):
        for f in job.output:
            self._record(self._metadata_path, {
                             "incomplete": True,
                             "external_jobid": external_jobid
                         }, f)

    def finished(self, job):
        version = str(
            job.rule.version) if job.rule.version is not None else None
        code = self._code(job.rule)
        input = self._input(job)
        log = self._log(job)
        params = self._params(job)
        shellcmd = job.shellcmd
        conda_env = self._conda_env(job)
        for f in job.expanded_output:
            rec_path = self._record_path(self._metadata_path, f)
            starttime = (os.path.getmtime(rec_path)
                         if os.path.exists(rec_path) else None)
            self._record(self._metadata_path, {
                "version": version,
                "code": code,
                "rule": job.rule.name,
                "input": input,
                "log": log,
                "params": params,
                "shellcmd": shellcmd,
                "incomplete": False,
                "starttime": starttime,
                "endtime": f.mtime,
                "job_hash": hash(job),
                "conda_env": conda_env,
                "singularity_img_url": job.singularity_img_url
            }, f)

    def cleanup(self, job):
        for f in job.expanded_output:
            self._delete_record(self._metadata_path, f)

    def incomplete(self, job):
        def marked_incomplete(f):
            return self._read_record(self._metadata_path, f).get("incomplete", False)

        return any(
            map(lambda f: f.exists and marked_incomplete(f), job.output))

    def external_jobids(self, job):
        return list(set(
            self._read_record(self._metadata_path, f)
                .get("external_jobid", None)
            for f in job.output))

    def metadata(self, path):
        return self._read_record(self._metadata_path, path)

    def version(self, path):
        return self.metadata(path).get("version")

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

    def version_changed(self, job, file=None):
        """Yields output files with changed versions of bool if file given."""
        return _bool_or_gen(self._version_changed, job, file=file)

    def code_changed(self, job, file=None):
        """Yields output files with changed code of bool if file given."""
        return _bool_or_gen(self._code_changed, job, file=file)

    def input_changed(self, job, file=None):
        """Yields output files with changed input of bool if file given."""
        return _bool_or_gen(self._input_changed, job, file=file)

    def params_changed(self, job, file=None):
        """Yields output files with changed params of bool if file given."""
        return _bool_or_gen(self._params_changed, job, file=file)

    def _version_changed(self, job, file=None):
        assert file is not None
        return self.version(file) != job.rule.version

    def _code_changed(self, job, file=None):
        assert file is not None
        return self.code(file) != self._code(job.rule)

    def _input_changed(self, job, file=None):
        assert file is not None
        return self.input(file) != self._input(job)

    def _params_changed(self, job, file=None):
        assert file is not None
        return self.params(file) != self._params(job)

    def noop(self, *args):
        pass

    def _b64id(self, s):
        return urlsafe_b64encode(str(s).encode()).decode()

    @lru_cache()
    def _code(self, rule):
        code = rule.run_func.__code__
        return b64encode(pickle_code(code)).decode()

    @lru_cache()
    def _conda_env(self, job):
        if job.conda_env:
            return b64encode(job.conda_env.content).decode()

    @lru_cache()
    def _input(self, job):
        return sorted(job.input)

    @lru_cache()
    def _log(self, job):
        return sorted(job.log)

    @lru_cache()
    def _params(self, job):
        return sorted(map(repr, job.params))

    @lru_cache()
    def _output(self, job):
        return sorted(job.output)

    def _record(self, subject, json_value, id):
        recpath = self._record_path(subject, id)
        os.makedirs(os.path.dirname(recpath), exist_ok=True)
        with open(recpath, "w") as f:
            json.dump(json_value, f)

    def _delete_record(self, subject, id):
        try:
            recpath = self._record_path(subject, id)
            os.remove(recpath)
            recdirs = os.path.relpath(os.path.dirname(recpath), start=subject)
            if recdirs != ".":
                os.removedirs(recdirs)
        except OSError as e:
            if e.errno != 2:  # not missing
                raise e

    @lru_cache()
    def _read_record_cached(self, subject, id):
        return self._read_record_uncached(subject, id)

    def _read_record_uncached(self, subject, id):
        if not self._exists_record(subject, id):
            return dict()
        with open(self._record_path(subject, id), "r") as f:
            return json.load(f)

    def _exists_record(self, subject, id):
        return os.path.exists(self._record_path(subject, id))

    def _locks(self, type):
        return (f for f, _ in listfiles(
            os.path.join(self._lockdir, "{{n,[0-9]+}}.{}.lock".format(type)))
                if not os.path.isdir(f))

    def _lock(self, files, type):
        for i in count(0):
            lockfile = os.path.join(self._lockdir,
                                    "{}.{}.lock".format(i, type))
            if not os.path.exists(lockfile):
                self._lockfile[type] = lockfile
                with open(lockfile, "w") as lock:
                    print(*files, sep="\n", file=lock)
                return

    def _record_path(self, subject, id):
        max_len = os.pathconf(
            subject,
            "PC_NAME_MAX") if os.name == "posix" else 255  # maximum NTFS and FAT32 filename length
        if max_len == 0:
            max_len = 255

        b64id = self._b64id(id)
        # split into chunks of proper length
        b64id = [b64id[i:i + max_len - 1]
                 for i in range(0, len(b64id), max_len - 1)]
        # prepend dirs with @ (does not occur in b64) to avoid conflict with b64-named files in the same dir
        b64id = ["@" + s for s in b64id[:-1]] + [b64id[-1]]
        path = os.path.join(subject, *b64id)
        return path

    def all_outputfiles(self):
        # we only look at output files that will be updated
        return jobfiles(self.dag.needrun_jobs, "output")

    def all_inputfiles(self):
        # we consider all input files, also of not running jobs
        return jobfiles(self.dag.jobs, "input")

    def deactivate_cache(self):
        self._read_record_cached.cache_clear()
        self._read_record = self._read_record_uncached


def _bool_or_gen(func, job, file=None):
    if file is None:
        return (f for f in job.expanded_output if func(job, file=f))
    else:
        return func(job, file=file)


def pickle_code(code):
    consts = [(pickle_code(const) if type(const) == type(code) else const)
              for const in code.co_consts]
    return pickle.dumps(
        (code.co_code, code.co_varnames, consts, code.co_names))
