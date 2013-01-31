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
        if not os.path.exists(self._started):
            os.mkdir(self._started)
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
        except OSError as e:
            if e.errno != 2:  # missing file
                raise e

    def started(self, job):
        with open(self.jobmarker(job), "w") as f:
            f.write("")

    def finished(self, job):
        try:
            os.remove(self.jobmarker(job))
        except OSError as e:
            if e.errno != 2:  # missing file
                raise e

    def incomplete(self, job):
        return os.access(self.jobmarker(job), os.W_OK)

    def noop(self):
        pass
    
    def jobmarker(self, job):
        return os.path.join(self._started, job.b64id)
