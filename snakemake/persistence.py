# -*- coding: utf-8 -*-

__author__ = "Johannes Köster"

import os
import shutil
import signal
import marshal
from base64 import urlsafe_b64encode
from functools import lru_cache, partial
from itertools import filterfalse, count

from snakemake.logging import logger
from snakemake.jobs import Job
from snakemake.utils import listfiles


class Persistence:

    def __init__(self, nolock=False, dag=None):
        self.path = os.path.abspath(".snakemake")
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        self._lockdir = os.path.join(self.path, "locks")
        if not os.path.exists(self._lockdir):
            os.mkdir(self._lockdir)

        self.dag = dag
        self._lockfile = dict()

        self._incomplete = os.path.join(self.path, "incomplete_files")
        self._version = os.path.join(self.path, "version_tracking")
        self._code = os.path.join(self.path, "code_tracking")
        self._rule = os.path.join(self.path, "rule_tracking")
        self._input = os.path.join(self.path, "input_tracking")
        self._params = os.path.join(self.path, "params_tracking")

        for d in (self._incomplete, self._version, self._code, self._rule, self._input, self._params):
            if not os.path.exists(d):
                os.mkdir(d)

        if nolock:
            self.lock = self.noop
            self.unlock = self.noop

    @property
    def files(self):
        if self._files is None:
            self._files = set(self.dag.output_files)
        return self._files

    @property
    def locked(self):
        strip = lambda lines: map(str.strip, lines)
        inputfiles = set(self.inputfiles())
        outputfiles = set(self.outputfiles())
        if os.path.exists(self._lockdir):
            for lockfile in self._locks("input"):
                with open(lockfile) as lock:
                    if not outputfiles.isdisjoint(strip(lock.readlines())):
                        return True
            for lockfile in self._locks("output"):
                with open(lockfile) as lock:
                    files = strip(lock.readlines())
                    if not outputfiles.isdisjoint(files):
                        return True
                    if not inputfiles.isdisjoint(files):
                        return True
        return False

    def lock(self):
        if self.locked:
            raise IOError("Another snakemake process "
                "has locked this directory.")
        self._lock(self.inputfiles(), "input")
        self._lock(self.outputfiles(), "output")

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
        self._delete_record(self._incomplete, path)
        self._delete_record(self._version, path)
        self._delete_record(self._code, path)
        self._delete_record(self._rule, path)
        self._delete_record(self._input, path)
        self._delete_record(self._params, path)

    def started(self, job):
        for f in job.output:
            self._record(self._incomplete, "", f)

    def finished(self, job):
        version = job.rule.version
        code = self.code(job.rule)
        input = self.input(job)
        params = self.params(job)
        for f in job.expanded_output:
            self._delete_record(self._incomplete, f)
            self._record(self._version, version, f)
            self._record(self._code, code, f, bin=True)
            self._record(self._rule, job.rule.name, f)
            self._record(self._input, input, f)
            self._record(self._params, params, f)

    def cleanup(self, job):
        for f in job.expanded_output:
            self._delete_record(self._incomplete, f)
            self._delete_record(self._version, f)
            self._delete_record(self._code, f)
            self._delete_record(self._rule, f)
            self._delete_record(self._input, f)
            self._delete_record(self._params, f)

    def incomplete(self, job):
        marked_incomplete = partial(self._exists_record, self._incomplete)
        return any(
            map(lambda f: f.exists and marked_incomplete(f), job.output))

    def version(self, path):
        if not os.path.exists(path):
            return None
        return self._read_record(self._version, path)

    def rule(self, path):
        if not os.path.exists(path):
            return None
        return self._read_record(self._rule, path)

    def version_changed(self, job, file=None):
        cr = partial(self._changed_records, self._version, job.rule.version)
        if file is None:
            return cr(*job.output)
        else:
            return bool(list(cr(file)))

    def code_changed(self, job, file=None):
        cr = partial(self._changed_records, self._code, self.code(job.rule), bin=True)
        if file is None:
            return cr(*job.output)
        else:
            return bool(list(cr(file)))

    def input_changed(self, job, file=None):
        cr = partial(self._changed_records, self._input, self.input(job))
        if file is None:
            return cr(*job.output)
        else:
            return bool(list(cr(file)))

    def params_changed(self, job, file=None):
        cr = partial(self._changed_records, self._params, self.params(job))
        if file is None:
            return cr(*job.output)
        else:
            return bool(list(cr(file)))

    def noop(self, *args):
        pass

    def b64id(self, s):
        return urlsafe_b64encode(str(s).encode()).decode()

    @lru_cache()
    def code(self, rule):
        return marshal.dumps(rule.run_func.__code__)

    @lru_cache()
    def input(self, job):
        return "\n".join(sorted(job.input))

    @lru_cache()
    def params(self, job):
        return "\n".join(sorted(job.params))

    @lru_cache()
    def output(self, job):
        return sorted(job.output)

    def _record(self, subject, value, id, bin=False):
        if value is not None:
            with open(
                os.path.join(subject, self.b64id(id)),
                "wb" if bin else "w") as f:
                f.write(value)

    def _delete_record(self, subject, id):
        try:
            os.remove(os.path.join(subject, self.b64id(id)))
        except OSError as e:
            if e.errno != 2:  # not missing
                raise e

    def _read_record(self, subject, id, bin=False):
        if not self._exists_record(subject, id):
            return None
        with open(
            os.path.join(subject, self.b64id(id)),
            "rb" if bin else "r") as f:
            return f.read()

    def _changed_records(self, subject, value, *ids, bin=False):
        equals = partial(self._equals_record, subject, value, bin=bin)
        return filter(
            lambda id: self._exists_record(subject, id) and not equals(id),
            ids)

    def _equals_record(self, subject, value, id, bin=False):
        return self._read_record(subject, id, bin=bin) == value

    def _exists_record(self, subject, id):
        return os.path.exists(os.path.join(subject, self.b64id(id)))

    def _locks(self, type):
        return (f for f, _ in listfiles(
            os.path.join(
                self._lockdir,
                "{{n,[0-9]+}}.{}.lock".format(type)))
            if not os.path.isdir(f))

    def _lock(self, files, type):
        for i in count(0):
            lockfile = os.path.join(
                self._lockdir, "{}.{}.lock".format(i, type))
            if not os.path.exists(lockfile):
                self._lockfile[type] = lockfile
                with open(lockfile, "w") as lock:
                    print(*files, sep="\n", file=lock)
                return

    def outputfiles(self):
        # we only look at output files that will be updated
        return Job.files(self.dag.needrun_jobs, "output")

    def inputfiles(self):
        # we consider all input files, also of not running jobs
        return Job.files(self.dag.jobs, "input")

