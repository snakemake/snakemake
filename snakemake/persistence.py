# -*- coding: utf-8 -*-

__author__ = "Johannes Köster"

import os
import signal


class Persistence:

    def __init__(self, nolock=False):
        self.path = os.path.abspath(".snakemake")
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        self._lock = os.path.join(self.path, "lock")
        self._started = os.path.join(self.path, "started_jobs")
        self._finished = os.path.join(self.path, "finished_jobs")
        if not os.path.exists(self._started):
            os.mkdir(self._started)
        if not os.path.exists(self._finished):
            os.mkdir(self._finished)
        if nolock:
            self.lock = self.noop
            self.unlock = self.noop

        for s in (signal.SIGTERM, signal.SIGABRT, signal.SIGINT):
            signal.signal(s, self.unlock)

    def lock(self):
        if os.path.exists(self._lock):
            raise IOError("Another snakemake process "
                "has locked this directory.")
        with open(self._lock, "w") as lock:
            lock.write("")

    def unlock(self):
        try:
            os.remove(self._lock)
        except OSError:
            pass

    def started(self, job):
        with open(os.path.join(self._started, job.b64id), "w") as f:
            f.write("")

    def finished(self, job):
        with open(os.path.join(self._finished, job.b64id), "w") as f:
            f.write("")

    def incomplete(self, job):
        return (os.path.exists(os.path.join(self._started, job.b64id))
            and not os.path.exists(os.path.join(self._finished, job.b64id)))

    def noop(self):
        pass
