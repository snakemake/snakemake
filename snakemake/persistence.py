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
        if not os.path.exists(self._incomplete):
            os.mkdir(self._incomplete)

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

    def started(self, job):
        for m in self.marker(job):
            self._mark_incomplete(m)

    def finished(self, job):
        for m in self.marker(job):
            self._mark_complete(m)

    def incomplete(self, job):
        return any(os.access(m, os.W_OK) for m in self.marker(job))

    def mark_complete(self, path):
        self._mark_complete(self.b64id(path))

    def noop(self):
        pass

    def marker(self, job):
        return map(self.b64id, job.output)

    @lru_cache()
    def b64id(self, path):
        return os.path.join(
            self._incomplete,
            urlsafe_b64encode(str(path).encode()).decode())

    def _mark_complete(self, marker):
        try:
            os.remove(marker)
        except OSError as e:
            if e.errno != 2:  # missing file
                raise e

    def _mark_incomplete(self, marker):
        with open(marker, "w") as m:
            m.write("")
