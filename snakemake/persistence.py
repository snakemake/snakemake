# -*- coding: utf-8 -*-

__author__ = "Johannes Köster"

import os
import signal
from base64 import urlsafe_b64encode
from functools import lru_cache


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

    def update_metadata(self, path):
        self._mark_complete(self.incomplete_marker(path))
        self.record_version(self.version_file(path))

    def started(self, job):
        for m in self.marker(job):
            self._mark_incomplete(m)

    def finished(self, job):
        for f in job.output:
            self.mark_complete(f)
            self.record_version(job, f)

    def incomplete(self, job):
        return any(os.access(m, os.W_OK)
            for m in map(self.incomplete_marker, job))

    def newversion(self, job):
        path = self.version_file(job)
        if os.path.exists(path):
            with open(path) as f:
                return job.rule.version != f.read()
        return False

    def record_version(self, job, path):
        if job.rule.version is not None:
            with open(self.version_file(path), "w") as f:
                f.write(str(job.rule.version))

    def noop(self):
        pass

    @lru_cache()
    def version_file(self, path):
        return os.path.join(self._version,
            self.b64id(path))

    @lru_cache()
    def incomplete_marker(self, path):
        return os.path.join(self._incomplete, self.b64id(path))

    def b64id(self, s):
        return urlsafe_b64encode(str(s).encode()).decode()

    def mark_complete(self, path):
        try:
            os.remove(self.incomplete_marker(path))
        except OSError as e:
            if e.errno != 2:  # missing file
                raise e

    def _mark_incomplete(self, marker):
        with open(marker, "w") as m:
            m.write("")
