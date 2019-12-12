__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2019, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import time
import csv
import json
from collections import defaultdict

import snakemake.jobs

fmt_time = time.ctime


class Stats:
    def __init__(self):
        self.starttime = dict()
        self.endtime = dict()

    def report_job_start(self, job):
        if job.is_group():
            for j in job:
                self.starttime[j] = time.time()
        else:
            self.starttime[job] = time.time()

    def report_job_end(self, job):
        if job.is_group():
            for j in job:
                self.endtime[j] = time.time()
        else:
            self.endtime[job] = time.time()

    @property
    def rule_stats(self):
        runtimes = defaultdict(list)
        for job, t in self.starttime.items():
            runtimes[job.rule].append(self.endtime[job] - t)
        for rule, runtimes in runtimes.items():
            yield (rule, sum(runtimes) / len(runtimes), min(runtimes), max(runtimes))

    @property
    def file_stats(self):
        for job, t in self.starttime.items():
            for f in job.expanded_output:
                start, stop = t, self.endtime[job]
                yield f, fmt_time(start), fmt_time(stop), stop - start, job

    @property
    def overall_runtime(self):
        if self.starttime and self.endtime:
            return max(self.endtime.values()) - min(self.starttime.values())
        else:
            return 0

    def to_json(self, path):
        rule_stats = {
            rule.name: {
                "mean-runtime": mean_runtime,
                "min-runtime": min_runtime,
                "max-runtime": max_runtime,
            }
            for rule, mean_runtime, min_runtime, max_runtime in self.rule_stats
        }
        file_stats = {
            f: {
                "start-time": start,
                "stop-time": stop,
                "duration": duration,
                "priority": job.priority
                if job.priority != snakemake.jobs.Job.HIGHEST_PRIORITY
                else "highest",
                "resources": dict(job.resources.items()),
            }
            for f, start, stop, duration, job in self.file_stats
        }

        with open(path, "w") as f:
            json.dump(
                {
                    "total_runtime": self.overall_runtime,
                    "rules": rule_stats,
                    "files": file_stats,
                },
                f,
                indent=4,
            )
