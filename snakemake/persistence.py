# -*- coding: utf-8 -*-

__author__ = "Johannes Köster"

import os
import signal
import marshal
import time
from base64 import urlsafe_b64encode
from functools import lru_cache, partial
from itertools import filterfalse


class Persistence:

    def __init__(self, nolock=False):
        self.path = os.path.abspath(".snakemake")
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        self._lock = os.path.join(self.path, "lock")

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
    def locked(self):
        return os.path.exists(self._lock)

    def lock(self):
        if self.locked:
            raise IOError("Another snakemake process "
                "has locked this directory.")
        with open(self._lock, "w") as lock:
            lock.write("")

    def unlock(self):
        try:
            os.remove(self._lock)
        except OSError as e:
            if e.errno != 2:  # missing file
                raise e

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

    def incomplete(self, job):
        return any(
            map(partial(self._exists_record, self._incomplete), job.output))

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
            if e.errno != 2: # not missing
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
