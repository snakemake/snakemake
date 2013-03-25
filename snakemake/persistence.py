# -*- coding: utf-8 -*-

__author__ = "Johannes Köster"

import os
import shutil
import signal
import marshal
from base64 import urlsafe_b64encode
from functools import lru_cache, partial
from itertools import filterfalse, count

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

        for d in (self._incomplete, self._version, self._code):
            if not os.path.exists(d):
                os.mkdir(d)

        if nolock:
            self.lock = self.noop
            self.unlock = self.noop

        for s in (signal.SIGTERM, signal.SIGABRT, signal.SIGINT):
            signal.signal(s, self.unlock)

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

    def unlock(self):
        for lockfile in self._lockfile.values():
            try:
                os.remove(lockfile)
            except OSError as e:
                if e.errno != 2:  # missing file
                    raise e

    def cleanup_locks(self):
        shutil.rmtree(self._lockdir)

    def cleanup_metadata(self, path):
        self._delete_record(self._incomplete, path)
        self._delete_record(self._version, path)
        self._delete_record(self._code, path)

    def started(self, job):
        for f in job.output:
            self._record(self._incomplete, "", f)

    def finished(self, job):
        version = job.rule.version
        code = self.code(job.rule)
        for f in job.output:
            self._delete_record(self._incomplete, f)
            self._record(self._version, version, f)
            self._record(self._code, code, f, bin=True)

    def cleanup(self, job):
        for f in job.output:
            self._delete_record(self._incomplete, f)
            self._delete_record(self._version, f)
            self._delete_record(self._code, f)

    def incomplete(self, job):
        marked_incomplete = partial(self._exists_record, self._incomplete)
        return any(
            map(lambda f: f.exists and marked_incomplete(f), job.output))

    def version(self, path):
        return self._read_record(self._version, path)

    def version_changed(self, job, file=None):
        if file is not None:
            return (self._exists_record(self._version, file)
                and not self._equals_record(
                    self._version, job.rule.version, file))
        return filterfalse(
            partial(self._equals_record, self._version, job.rule.version),
            job.output)

    def code_changed(self, job, file=None):
        if file is not None:
            return (self._exists_record(self._code, file)
                and not self._equals_record(
                    self._code, self.code(job.rule), file, bin=True))
        return filterfalse(
            partial(
                self._equals_record, self._code,
                self.code(job.rule), bin=True),
            job.output)

    def noop(self):
        pass

    def b64id(self, s):
        return urlsafe_b64encode(str(s).encode()).decode()

    @lru_cache()
    def code(self, rule):
        return marshal.dumps(rule.run_func.__code__)

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

