__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import shutil
import signal
import marshal
import pickle
from base64 import urlsafe_b64encode
from functools import lru_cache, partial
from itertools import filterfalse, count

from snakemake.logging import logger
from snakemake.jobs import jobfiles
from snakemake.utils import listfiles


class Persistence:
    def __init__(self, nolock=False, dag=None, warn_only=False):
        self.path = os.path.abspath(".snakemake")
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        self._lockdir = os.path.join(self.path, "locks")
        if not os.path.exists(self._lockdir):
            os.mkdir(self._lockdir)

        self.dag = dag
        self._lockfile = dict()

        self._incomplete_path = os.path.join(self.path, "incomplete_files")
        self._version_path = os.path.join(self.path, "version_tracking")
        self._code_path = os.path.join(self.path, "code_tracking")
        self._rule_path = os.path.join(self.path, "rule_tracking")
        self._input_path = os.path.join(self.path, "input_tracking")
        self._params_path = os.path.join(self.path, "params_tracking")
        self._shellcmd_path = os.path.join(self.path, "shellcmd_tracking")
        self.shadow_path = os.path.join(self.path, "shadow")

        for d in (self._incomplete_path, self._version_path, self._code_path,
                  self._rule_path, self._input_path, self._params_path,
                  self._shellcmd_path, self.shadow_path):
            if not os.path.exists(d):
                os.mkdir(d)

        if nolock:
            self.lock = self.noop
            self.unlock = self.noop
        if warn_only:
            self.lock = self.lock_warn_only
            self.unlock = self.noop

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
                "means that another Snakemake instance is running on this directory."
                "Another possiblity is that a previous run exited unexpectedly.")

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
        self._delete_record(self._incomplete_path, path)
        self._delete_record(self._version_path, path)
        self._delete_record(self._code_path, path)
        self._delete_record(self._rule_path, path)
        self._delete_record(self._input_path, path)
        self._delete_record(self._params_path, path)
        self._delete_record(self._shellcmd_path, path)

    def cleanup_shadow(self):
        if os.path.exists(self.shadow_path):
            shutil.rmtree(self.shadow_path)
            os.mkdir(self.shadow_path)

    def started(self, job):
        for f in job.output:
            self._record(self._incomplete_path, "", f)

    def finished(self, job):
        version = str(
            job.rule.version) if job.rule.version is not None else None
        code = self._code(job.rule)
        input = self._input(job)
        params = self._params(job)
        shellcmd = self._shellcmd(job)
        for f in job.expanded_output:
            self._delete_record(self._incomplete_path, f)
            self._record(self._version_path, version, f)
            self._record(self._code_path, code, f, bin=True)
            self._record(self._rule_path, job.rule.name, f)
            self._record(self._input_path, input, f)
            self._record(self._params_path, params, f)
            self._record(self._shellcmd_path, shellcmd, f)

    def cleanup(self, job):
        for f in job.expanded_output:
            self._delete_record(self._incomplete_path, f)
            self._delete_record(self._version_path, f)
            self._delete_record(self._code_path, f)
            self._delete_record(self._rule_path, f)
            self._delete_record(self._input_path, f)
            self._delete_record(self._params_path, f)
            self._delete_record(self._shellcmd_path, f)

    def incomplete(self, job):
        marked_incomplete = partial(self._exists_record, self._incomplete_path)
        return any(
            map(lambda f: f.exists and marked_incomplete(f), job.output))

    def version(self, path):
        return self._read_record(self._version_path, path)

    def rule(self, path):
        return self._read_record(self._rule_path, path)

    def input(self, path):
        files = self._read_record(self._input_path, path)
        if files is not None:
            return files.split("\n")
        return None

    def shellcmd(self, path):
        return self._read_record(self._shellcmd_path, path)

    def version_changed(self, job, file=None):
        cr = partial(self._changed_records, self._version_path,
                     job.rule.version)
        if file is None:
            return cr(*job.output)
        else:
            return bool(list(cr(file)))

    def code_changed(self, job, file=None):
        cr = partial(self._changed_records, self._code_path,
                     self._code(job.rule),
                     bin=True)
        if file is None:
            return cr(*job.output)
        else:
            return bool(list(cr(file)))

    def input_changed(self, job, file=None):
        cr = partial(self._changed_records, self._input_path, self._input(job))
        if file is None:
            return cr(*job.output)
        else:
            return bool(list(cr(file)))

    def params_changed(self, job, file=None):
        cr = partial(self._changed_records, self._params_path,
                     self._params(job))
        if file is None:
            return cr(*job.output)
        else:
            return bool(list(cr(file)))

    def noop(self, *args):
        pass

    def _b64id(self, s):
        return urlsafe_b64encode(str(s).encode()).decode()

    @lru_cache()
    def _code(self, rule):
        code = rule.run_func.__code__
        return pickle_code(code)

    @lru_cache()
    def _input(self, job):
        return "\n".join(sorted(job.input))

    @lru_cache()
    def _params(self, job):
        return "\n".join(sorted(map(repr, job.params)))

    @lru_cache()
    def _output(self, job):
        return sorted(job.output)

    @lru_cache()
    def _shellcmd(self, job):
        return job.shellcmd

    def _record(self, subject, value, id, bin=False):
        recpath = self._record_path(subject, id)
        if value is not None:
            os.makedirs(os.path.dirname(recpath), exist_ok=True)
            with open(recpath, "wb" if bin else "w") as f:
                f.write(value)
        else:
            if os.path.exists(recpath):
                os.remove(recpath)

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

    def _read_record(self, subject, id, bin=False):
        if not self._exists_record(subject, id):
            return None
        with open(self._record_path(subject, id), "rb" if bin else "r") as f:
            return f.read()

    def _changed_records(self, subject, value, *ids, bin=False):
        equals = partial(self._equals_record, subject, value, bin=bin)
        return filter(
            lambda id: self._exists_record(subject, id) and not equals(id),
            ids)

    def _equals_record(self, subject, value, id, bin=False):
        return self._read_record(subject, id, bin=bin) == value

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


def pickle_code(code):
    consts = [(pickle_code(const) if type(const) == type(code) else const)
              for const in code.co_consts]
    return pickle.dumps(
        (code.co_code, code.co_varnames, consts, code.co_names))
