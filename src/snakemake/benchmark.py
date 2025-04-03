__author__ = "Manuel Holtgrewe"
__copyright__ = "Copyright 2022, Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"
__license__ = "MIT"

import contextlib
import datetime
from itertools import chain
import os
import time
import threading
from pathlib import Path

from snakemake.logging import logger

#: Interval (in seconds) between measuring resource usage
BENCHMARK_INTERVAL = 30
#: Interval (in seconds) between measuring resource usage before
#: BENCHMARK_INTERVAL
BENCHMARK_INTERVAL_SHORT = 0.5


class BenchmarkRecord:
    """Record type for benchmark times"""

    @classmethod
    def get_header(klass, extended_fmt=False):
        header = [
            "s",
            "h:m:s",
            "max_rss",
            "max_vms",
            "max_uss",
            "max_pss",
            "io_in",
            "io_out",
            "mean_load",
            "cpu_time",
        ]

        if extended_fmt:
            header += [
                "jobid",
                "rule_name",
                "wildcards",
                "params",
                "threads",
                "cpu_usage",
                "resources",
                "input_size_mb",
            ]

        return header

    def __init__(
        self,
        jobid=None,
        rule_name=None,
        wildcards=None,
        params=None,
        running_time=None,
        max_rss=None,
        max_vms=None,
        max_uss=None,
        max_pss=None,
        io_in=None,
        io_out=None,
        cpu_usage=None,
        cpu_time=None,
        resources=None,
        threads=None,
        input=None,
    ):
        #: Job ID
        self.jobid = (jobid,)
        #: Rule name
        self.rule_name = (rule_name,)
        #: Job wildcards
        self.wildcards = (wildcards,)
        #: Job parameters
        self.params = (params,)
        #: Running time in seconds
        self.running_time = running_time
        #: Maximal RSS in MB
        self.max_rss = max_rss
        #: Maximal VMS in MB
        self.max_vms = max_vms
        #: Maximal USS in MB
        self.max_uss = max_uss
        #: Maximal PSS in MB
        self.max_pss = max_pss
        #: I/O read in bytes
        self.io_in = io_in
        #: I/O written in bytes
        self.io_out = io_out
        #: Count of CPU seconds, divide by running time to get mean load estimate
        self.cpu_usage = cpu_usage or 0
        #: CPU usage (user and system) in seconds
        self.cpu_time = cpu_time or 0
        #: Job resources
        self.resources = (resources,)
        #: Job threads
        self.threads = (threads,)
        #: Job input
        self.input = input
        #: First time when we measured CPU load, for estimating total running time
        self.first_time = None
        #: Previous point when measured CPU load, for estimating total running time
        self.prev_time = None
        #: Set with procs that has been skipped
        self.processed_procs = dict()
        #: Set with procs that has been saved
        self.skipped_procs = set()
        #: Track if data has been collected
        self.data_collected = False

    def timedelta_to_str(self, x):
        """Conversion of timedelta to str without fractions of seconds"""
        mm, ss = divmod(x.seconds, 60)
        hh, mm = divmod(mm, 60)
        s = "%d:%02d:%02d" % (hh, mm, ss)
        if x.days:

            def plural(n):
                return n, abs(n) != 1 and "s" or ""

            s = ("%d day%s, " % plural(x.days)) + s
        return s

    def mean_load(self):
        return self.cpu_usage / self.running_time

    def parse_wildcards(self):
        return {key: value for key, value in self.wildcards.items()}

    def parse_params(self):
        return {key: value for key, value in self.params.items()}

    def parse_resources(self):
        return {key: value for key, value in self.resources.items()}

    def input_size_mb(self):
        return {file: Path(file).stat().st_size / 1024 / 1024 for file in self.input}

    def get_benchmarks(self, extended_fmt=False):
        logger.debug(
            f"Stats included in benchmarks file: {self.get_header(extended_fmt)}"
        )
        if self.skipped_procs:
            logger.debug(
                "Benchmark: not collected for "
                "; ".join(
                    [
                        f"{{'pid': {record[0]}, 'name': '{record[1]}''}}"
                        for record in self.skipped_procs
                    ]
                )
            )
            logger.debug(
                "Benchmark: collected for "
                "; ".join(
                    [
                        f"{{'pid': {record[0]}, 'name': '{record[1]}'}}"
                        for record in self.processed_procs
                    ]
                )
            )

        # If no data has been collect mem and cpu statistics will be printed as NA
        # to make it possible to distinguish this case from processes that complete instantly
        if not self.data_collected:
            logger.warning(
                "Benchmark: unable to collect cpu and memory benchmark statistics"
            )
        record = [
            f"{self.running_time:.4f}",
            self.timedelta_to_str(datetime.timedelta(seconds=self.running_time)),
            self.max_rss if self.data_collected else "NA",
            self.max_vms if self.data_collected else "NA",
            self.max_uss if self.data_collected else "NA",
            self.max_pss if self.data_collected else "NA",
            self.io_in if self.data_collected else "NA",
            self.io_out if self.data_collected else "NA",
            self.mean_load() if self.data_collected else "NA",
            self.cpu_time if self.data_collected else "NA",
        ]
        if extended_fmt:
            record += [
                self.jobid,
                self.rule_name,
                self.parse_wildcards(),
                self.parse_params(),
                self.threads,
                self.cpu_usage if self.data_collected else "NA",
                self.parse_resources(),
                self.input_size_mb(),
            ]
        return record

    def to_tsv(self, extended_fmt):
        """Return ``str`` with the TSV representation of this record"""

        def to_tsv_str(x):
            """Conversion of value to str for TSV (None becomes "-")"""
            if x is None:
                return "-"
            elif isinstance(x, float):
                return f"{x:.2f}"
            else:
                return str(x)

        return "\t".join(map(to_tsv_str, self.get_benchmarks(extended_fmt)))

    def to_json(self, extended_fmt):
        """Return ``str`` with the JSON representation of this record"""
        import json

        return json.dumps(
            dict(zip(self.get_header(extended_fmt), self.get_benchmarks(extended_fmt))),
            sort_keys=True,
        )


class DaemonTimer(threading.Thread):
    """A variant of threading.The timer that is daemonized"""

    def __init__(self, interval, function, args=None, kwargs=None):
        threading.Thread.__init__(self, daemon=True)
        self.interval = interval
        self.function = function
        self.args = args if args is not None else []
        self.kwargs = kwargs if kwargs is not None else {}
        self.finished = threading.Event()

    def cancel(self):
        """Stop the timer if it hasn't finished yet."""
        self.finished.set()

    def run(self):
        self.finished.wait(self.interval)
        if not self.finished.is_set():
            self.function(*self.args, **self.kwargs)
        self.finished.set()


class ScheduledPeriodicTimer:
    """Scheduling of periodic events

    Up to self._interval, schedule actions per second, above schedule events
    in self._interval second gaps.
    """

    def __init__(self, interval):
        self._times_called = 0
        self._interval = interval
        self._timer = None
        self._stopped = True

    def start(self):
        """Start the intervalic timer"""
        self.work()
        self._times_called += 1
        self._stopped = False
        if self._times_called > self._interval:
            self._timer = DaemonTimer(self._interval, self._action)
        else:
            self._timer = DaemonTimer(BENCHMARK_INTERVAL_SHORT, self._action)
        self._timer.start()

    def _action(self):
        """Internally, called by timer"""
        self.work()
        self._times_called += 1
        if self._times_called > self._interval:
            self._timer = DaemonTimer(self._interval, self._action)
        else:
            self._timer = DaemonTimer(BENCHMARK_INTERVAL_SHORT, self._action)
        self._timer.start()

    def work(self):
        """Override to perform the action"""
        raise NotImplementedError("Override me!")

    def cancel(self):
        """Call to cancel any events"""
        self._timer.cancel()
        self._stopped = True


class BenchmarkTimer(ScheduledPeriodicTimer):
    """Allows easy observation of a given PID for resource usage"""

    def __init__(self, pid, bench_record, interval=BENCHMARK_INTERVAL):
        import psutil

        ScheduledPeriodicTimer.__init__(self, interval)
        #: PID of observed process
        self.pid = pid
        self.main = psutil.Process(self.pid)
        #: ``BenchmarkRecord`` to write results to
        self.bench_record = bench_record
        #: Cache of processes to keep track of cpu percent
        self.procs = {}

    def work(self):
        """Write statistics"""
        import psutil

        try:
            self._update_record()
        except psutil.NoSuchProcess:
            pass  # skip, process died in flight
        except AttributeError:
            pass  # skip, process died in flight

    def _update_record(self):
        """Perform the actual measurement"""
        import psutil

        # Memory measurements
        rss, vms, uss, pss = 0, 0, 0, 0
        # I/O measurements
        io_in, io_out = 0, 0
        check_io = True
        # CPU seconds
        cpu_usage = 0
        # CPU usage time
        cpu_time = 0

        data_collected = False
        # Iterate over process and all children
        try:
            this_time = time.time()
            for proc in chain((self.main,), self.main.children(recursive=True)):
                proc = self.procs.setdefault(proc.pid, proc)
                with proc.oneshot():
                    if self.bench_record.prev_time:
                        cpu_usage += proc.cpu_percent() * (
                            this_time - self.bench_record.prev_time
                        )
                    # Makes it possible to summarize information about the process even
                    # if the benchmark has tried to access a process that the user does
                    # not have access to.
                    try:
                        meminfo = proc.memory_full_info()
                    except psutil.Error:
                        # Continue to fetch information about the remaining processes
                        # save skipped processes pid and name for debugging
                        self.bench_record.skipped_procs.add((proc.pid, proc.name()))
                        continue
                    rss += meminfo.rss
                    vms += meminfo.vms
                    uss += meminfo.uss
                    pss += meminfo.pss

                    if check_io:
                        try:
                            ioinfo = proc.io_counters()
                            io_in += ioinfo.read_bytes
                            io_out += ioinfo.write_bytes
                        except NotImplementedError as nie:
                            # OS doesn't track IO
                            check_io = False

                    cpu_times = proc.cpu_times()
                    self.bench_record.processed_procs[(proc.pid, proc.name())] = (
                        cpu_times.user + cpu_times.system
                    )

            cpu_time = sum(self.bench_record.processed_procs.values())

            self.bench_record.prev_time = this_time
            if not self.bench_record.first_time:
                self.bench_record.prev_time = this_time

            rss /= 1024 * 1024
            vms /= 1024 * 1024
            uss /= 1024 * 1024
            pss /= 1024 * 1024

            if check_io:
                io_in /= 1024 * 1024
                io_out /= 1024 * 1024
            else:
                io_in = None
                io_out = None
            data_collected = True
        except psutil.Error as e:
            return

        # Update benchmark record's RSS and VMS
        if data_collected:
            self.bench_record.data_collected = True
            self.bench_record.max_rss = max(self.bench_record.max_rss or 0, rss)
            self.bench_record.max_vms = max(self.bench_record.max_vms or 0, vms)
            self.bench_record.max_uss = max(self.bench_record.max_uss or 0, uss)
            self.bench_record.max_pss = max(self.bench_record.max_pss or 0, pss)

            self.bench_record.io_in = io_in
            self.bench_record.io_out = io_out

            self.bench_record.cpu_usage += cpu_usage
            self.bench_record.cpu_time = cpu_time


@contextlib.contextmanager
def benchmarked(pid=None, benchmark_record=None, interval=BENCHMARK_INTERVAL):
    """Measure benchmark parameters while within the context manager

    Yields a ``BenchmarkRecord`` with the results (values are set after
    leaving context).

    If ``pid`` is ``None`` then the PID of the current process will be used.
    If ``benchmark_record`` is ``None`` then a new ``BenchmarkRecord`` is
    created and returned, otherwise, the object passed as this parameter is
    returned.

    Usage::

        with benchmarked() as bench_result:
            pass
    """
    result = benchmark_record or BenchmarkRecord()
    if pid is False:
        yield result
    else:
        start_time = time.time()
        bench_thread = BenchmarkTimer(int(pid or os.getpid()), result, interval)
        bench_thread.start()
        yield result
        bench_thread.cancel()
        result.running_time = time.time() - start_time


def print_benchmark_tsv(records, file_, extended_fmt):
    """Write benchmark records to file-like the object"""
    logger.debug("Benchmarks in TSV format")
    print("\t".join(BenchmarkRecord.get_header(extended_fmt)), file=file_)
    for r in records:
        print(r.to_tsv(extended_fmt), file=file_)


def print_benchmark_jsonl(records, file_, extended_fmt):
    """Write benchmark records to file-like the object"""
    logger.debug("Benchmarks in JSONL format")
    for r in records:
        print(r.to_json(extended_fmt), file=file_)


def write_benchmark_records(records, path, extended_fmt):
    """Write benchmark records to file at path"""
    with open(path, "wt") as f:
        if path.endswith(".jsonl"):
            print_benchmark_jsonl(records, f, extended_fmt)
        else:
            print_benchmark_tsv(records, f, extended_fmt)
