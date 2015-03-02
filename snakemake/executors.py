
import os
import time
import datetime
import json
import textwrap
import stat
import shutil
import random
import string
import threading
import concurrent.futures
import subprocess
import signal
from functools import partial
from itertools import chain

from snakemake.jobs import Job
from snakemake.shell import shell
from snakemake.logging import logger
from snakemake.stats import Stats
from snakemake.utils import format, Unformattable
from snakemake.io import get_wildcard_names
from snakemake.exceptions import print_exception, get_exception_origin
from snakemake.exceptions import format_error, RuleException
from snakemake.exceptions import ClusterJobException, ProtectedOutputException, WorkflowError
from snakemake.futures import ProcessPoolExecutor


class AbstractExecutor:

    def __init__(self, workflow, dag,
        printreason=False, quiet=False,
        printshellcmds=False, printthreads=True, latency_wait=3,
        benchmark_repeats=3):
        self.workflow = workflow
        self.dag = dag
        self.quiet = quiet
        self.printreason = printreason
        self.printshellcmds = printshellcmds
        self.printthreads = printthreads
        self.latency_wait = latency_wait
        self.benchmark_repeats = benchmark_repeats

    def run(
        self, job, callback=None, submit_callback=None, error_callback=None):
        job.check_protected_output()
        self._run(job)
        callback(job)

    def shutdown(self):
        pass

    def _run(self, job):
        self.printjob(job)

    def rule_prefix(self, job):
        return "local " if self.workflow.is_local(job.rule) else ""

    def printjob(self, job):
        # skip dynamic jobs that will be "executed" only in dryrun mode
        if self.dag.dynamic(job):
            return

        def format_files(job, io, ruleio, dynamicio):
            for f in io:
                f_ = ruleio[f]
                if f in dynamicio:
                    yield "{} (dynamic)".format(f.format_dynamic())
                else:
                    yield f

        priority = self.dag.priority(job)
        logger.job_info(
            jobid=self.dag.jobid(job),
            msg=job.message,
            name=job.rule.name,
            local=self.workflow.is_local(job.rule),
            input=list(format_files(job, job.input, job.ruleio, job.dynamic_input)),
            output=list(format_files(job, job.output, job.ruleio, job.dynamic_output)),
            log=job.log,
            benchmark=job.benchmark,
            reason=str(self.dag.reason(job)),
            resources=job.resources_dict,
            priority="highest" if priority == Job.HIGHEST_PRIORITY else priority,
            threads=job.threads)

        if job.dynamic_output:
            logger.info(
                    "Subsequent jobs will be added dynamically "
                    "depending on the output of this rule")

    def print_job_error(self, job):
        logger.error("Error in job {} while creating output file{} {}.".format(job, "s" if len(job.output) > 1 else "", ", ".join(job.output)))

    def finish_job(self, job):
        self.dag.handle_touch(job)
        self.dag.check_output(job, wait=self.latency_wait)
        self.dag.handle_protected(job)
        self.dag.handle_temp(job)


class DryrunExecutor(AbstractExecutor):

    def _run(self, job):
        super()._run(job)
        logger.shellcmd(job.shellcmd)


class RealExecutor(AbstractExecutor):

    def __init__(
        self, workflow, dag,
        printreason=False, quiet=False, printshellcmds=False, latency_wait=3,
        benchmark_repeats=3):
        super().__init__(
            workflow, dag, printreason=printreason,
            quiet=quiet, printshellcmds=printshellcmds,
            latency_wait=latency_wait,
            benchmark_repeats=benchmark_repeats)
        self.stats = Stats()

    def _run(self, job, callback=None, error_callback=None):
        super()._run(job)
        self.stats.report_job_start(job)
        try:
            self.workflow.persistence.started(job)
        except IOError as e:
            logger.info("Failed to set marker file for job started ({}). "
                "Snakemake will work, but cannot ensure that output files "
                "are complete in case of a kill signal or power loss. "
                "Please ensure write permissions for the "
                "directory {}".format(
                    e, self.workflow.persistence.path))

    def finish_job(self, job):
        super().finish_job(job)
        self.stats.report_job_end(job)
        try:
            self.workflow.persistence.finished(job)
        except IOError as e:
            logger.info("Failed to remove marker file for job started "
                "({}). Please ensure write permissions for the "
                "directory {}".format(
                    e, self.workflow.persistence.path))


class TouchExecutor(RealExecutor):

    def run(
        self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)
        try:
            for f in job.expanded_output:
                f.touch()
            if job.benchmark:
                job.benchmark.touch()
            time.sleep(0.1)
            self.finish_job(job)
            callback(job)
        except OSError as ex:
            print_exception(ex, self.workflow.linemaps)
            error_callback(job)


_ProcessPoolExceptions = (KeyboardInterrupt,)
try:
    from concurrent.futures.process import BrokenProcessPool
    _ProcessPoolExceptions = (KeyboardInterrupt, BrokenProcessPool)
except ImportError:
    pass


class CPUExecutor(RealExecutor):

    def __init__(
        self,
        workflow,
        dag,
        cores,
        printreason=False,
        quiet=False,
        printshellcmds=False,
        threads=False,
        latency_wait=3,
        benchmark_repeats=3
    ):
        super().__init__(
            workflow, dag, printreason=printreason, quiet=quiet,
            printshellcmds=printshellcmds, latency_wait=latency_wait,
            benchmark_repeats=benchmark_repeats)

        self.pool = (concurrent.futures.ThreadPoolExecutor(max_workers=cores)
            if threads
            else ProcessPoolExecutor(max_workers=cores))

    def run(
        self, job, callback=None, submit_callback=None, error_callback=None):
        job.prepare()
        super()._run(job)

        log = None
        if job.log is not None:
            log = str(job.log)
        benchmark = None
        if job.benchmark is not None:
            benchmark = str(job.benchmark)

        future = self.pool.submit(
            run_wrapper, job.rule.run_func, job.input.plainstrings(), job.output.plainstrings(), job.params,
            job.wildcards, job.threads, job.resources, log, job.rule.version, benchmark, self.benchmark_repeats, self.workflow.linemaps)
        future.add_done_callback(partial(
            self._callback, job, callback, error_callback))

    def shutdown(self):
        self.pool.shutdown()

    def cancel(self):
        self.pool.shutdown()

    def _callback(self, job, callback, error_callback, future):
        try:
            ex = future.exception()
            if ex:
                raise ex
            self.finish_job(job)
            callback(job)
        except _ProcessPoolExceptions:
            job.cleanup()
            self.workflow.persistence.cleanup(job)
            # no error callback, just silently ignore the interrupt as the main scheduler is also killed
        except (Exception, BaseException) as ex:
            self.print_job_error(job)
            print_exception(ex, self.workflow.linemaps)
            job.cleanup()
            self.workflow.persistence.cleanup(job)
            error_callback(job)


class ClusterExecutor(RealExecutor):

    default_jobscript = "jobscript.sh"

    def __init__(
        self,
        workflow,
        dag,
        cores,
        jobname="snakejob.{rulename}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        latency_wait=3,
        benchmark_repeats=3
    ):
        super().__init__(
            workflow, dag, printreason=printreason, quiet=quiet,
            printshellcmds=printshellcmds, latency_wait=latency_wait)
        if workflow.snakemakepath is None:
            raise ValueError(
            "Cluster executor needs to know the path "
            "to the snakemake binary.")

        jobscript = workflow.jobscript
        if jobscript is None:
            jobscript = os.path.join(
                os.path.dirname(__file__),
                self.default_jobscript)
        try:
            with open(jobscript) as f:
                self.jobscript = f.read()
        except IOError as e:
            raise WorkflowError(e)

        if not "jobid" in get_wildcard_names(jobname):
            raise WorkflowError(
                "Defined jobname (\"{}\") has to contain the wildcard {jobid}.")

        self.exec_job = (
            'cd {workflow.workdir_init} && '
            '{workflow.snakemakepath} --snakefile {workflow.snakefile} '
            '--force -j{cores} --keep-target-files '
            '--wait-for-files {job.input} --latency-wait {latency_wait} '
            '--benchmark-repeats {benchmark_repeats} '
            '{overwrite_workdir} {overwrite_config} --nocolor '
            '--notemp --quiet --nolock {target}'
        )

        if printshellcmds:
            self.exec_job += " --printshellcmds "

        if not any(dag.dynamic_output_jobs):
            # disable restiction to target rule in case of dynamic rules!
            self.exec_job += " --allowed-rules {job.rule.name} "
        self.jobname = jobname
        self.threads = []
        self._tmpdir = None
        self.cores = cores if cores else ""

    def shutdown(self):
        for thread in self.threads:
            thread.join()
        shutil.rmtree(self.tmpdir)

    def cancel(self):
        self.shutdown()    
    
    @property
    def tmpdir(self):
        if self._tmpdir is None:
            while True:
                self._tmpdir = ".snakemake/tmp." + "".join(random.sample(
                    string.ascii_uppercase + string.digits, 6))
                if not os.path.exists(self._tmpdir):
                    os.mkdir(self._tmpdir)
                    break
        return os.path.abspath(self._tmpdir)

    def get_jobscript(self, job):
        return os.path.join(self.tmpdir, self.jobname.format(rulename=job.rule.name, jobid=self.dag.jobid(job)))

    def spawn_jobscript(self, job, jobscript, **kwargs):
        overwrite_workdir = ""
        if self.workflow.overwrite_workdir:
            overwrite_workdir = "--directory {}".format(self.workflow.overwrite_workdir)
        overwrite_config = ""
        if self.workflow.overwrite_configfile:
            overwrite_config = "--configfile {}".format(self.workflow.overwrite_configfile)
        if self.workflow.config_args:
            overwrite_config += "--config {}".format(" ".join(self.workflow.config_args))

        target = job.output if job.output else job.rule.name
        format = partial(
            str.format,
            job=job,
            overwrite_workdir=overwrite_workdir,
            overwrite_config=overwrite_config,
            workflow=self.workflow,
            cores=self.cores,
            properties=job.json(),
            latency_wait=self.latency_wait,
            benchmark_repeats=self.benchmark_repeats,
            target=target,
            **kwargs)
        try:
            exec_job = format(self.exec_job)
            with open(jobscript, "w") as f:
                print(format(self.jobscript, exec_job=exec_job), file=f)
        except KeyError as e:
            raise WorkflowError(
                "Error formatting jobscript: {} not found\n"
                "Make sure that your custom jobscript it up to date.".format(e))
        os.chmod(jobscript, os.stat(jobscript).st_mode | stat.S_IXUSR)


class GenericClusterExecutor(ClusterExecutor):

    def __init__(
        self,
        workflow,
        dag,
        cores,
        submitcmd="qsub",
        jobname="snakejob.{rulename}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        latency_wait=3,
        benchmark_repeats=3
    ):
        super().__init__(
            workflow, dag, cores, jobname=jobname,
            printreason=printreason, quiet=quiet,
            printshellcmds=printshellcmds, latency_wait=latency_wait,
            benchmark_repeats=benchmark_repeats)
        self.submitcmd = submitcmd
        self.external_jobid = dict()
        self.exec_job += ' && touch "{jobfinished}" || touch "{jobfailed}"'

    def cancel(self):
        logger.info("Will exit after finishing currently running jobs.")
        self.shutdown()

    def run(
        self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)
        workdir = os.getcwd()
        jobid = self.dag.jobid(job)

        jobscript = self.get_jobscript(job)
        jobfinished = os.path.join(self.tmpdir, "{}.jobfinished".format(jobid))
        jobfailed = os.path.join(self.tmpdir, "{}.jobfailed".format(jobid))
        self.spawn_jobscript(
            job, jobscript,
            jobfinished=jobfinished, jobfailed=jobfailed)

        deps = " ".join(self.external_jobid[f] for f in job.input if f in self.external_jobid)
        submitcmd = job.format_wildcards(self.submitcmd, dependencies=deps)
        try:
            ext_jobid = subprocess.check_output(
                '{submitcmd} "{jobscript}"'.format(
                    submitcmd=submitcmd,
                    jobscript=jobscript),
                shell=True).decode().split("\n")
        except subprocess.CalledProcessError as ex:
            raise WorkflowError("Error executing jobscript (exit code {}):\n{}".format(ex.returncode, ex.output.decode()), rule=job.rule)
        if ext_jobid and ext_jobid[0]:
            ext_jobid = ext_jobid[0]
            self.external_jobid.update((f, ext_jobid) for f in job.output)
            logger.debug("Submitted job {} with external jobid {}.".format(jobid, ext_jobid))

        thread = threading.Thread(
            target=self._wait_for_job,
            args=(job, callback, error_callback,
                jobscript, jobfinished, jobfailed))
        thread.daemon = True
        thread.start()
        self.threads.append(thread)

        submit_callback(job)

    def _wait_for_job(
        self, job, callback, error_callback,
        jobscript, jobfinished, jobfailed):
        while True:
            if os.path.exists(jobfinished):
                os.remove(jobfinished)
                os.remove(jobscript)
                self.finish_job(job)
                callback(job)
                return
            if os.path.exists(jobfailed):
                os.remove(jobfailed)
                os.remove(jobscript)
                self.print_job_error(job)
                print_exception(
                    ClusterJobException(job, self.dag.jobid(job), self.get_jobscript(job)), self.workflow.linemaps)
                error_callback(job)
                return
            time.sleep(1)


class DRMAAExecutor(ClusterExecutor):

    def __init__(
        self,
        workflow,
        dag,
        cores,
        jobname="snakejob.{rulename}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        drmaa_args="",
        latency_wait=3,
        benchmark_repeats=3
    ):
        super().__init__(workflow, dag, cores, jobname=jobname,
            printreason=printreason, quiet=quiet,
            printshellcmds=printshellcmds, latency_wait=latency_wait,
            benchmark_repeats=benchmark_repeats)
        try:
            import drmaa
        except ImportError:
            raise WorkflowError("Python support for DRMAA is not installed. "
            "Please install it, e.g. with easy_install3 --user drmaa")
        except RuntimeError as e:
            raise WorkflowError("Error loading drmaa support:\n{}".format(e))
        self.session = drmaa.Session()
        self.drmaa_args=drmaa_args
        self.session.initialize()
        self.submitted = list()

    def cancel(self):
        from drmaa.const import JobControlAction
        for jobid in self.submitted:
            self.session.control(jobid, JobControlAction.TERMINATE)
        self.shutdown()

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)
        jobscript = self.get_jobscript(job)
        self.spawn_jobscript(job, jobscript)
        import drmaa
        try:
            jt = self.session.createJobTemplate()
            jt.remoteCommand = jobscript
            jt.nativeSpecification = job.format_wildcards(self.drmaa_args)

            jobid = self.session.runJob(jt)
        except (drmaa.errors.InternalException, drmaa.errors.InvalidAttributeValueException) as e:
            print_exception(WorkflowError("DRMAA Error: {}".format(e)), self.workflow.linemaps)
            error_callback(job)
            return
        logger.info("Submitted DRMAA job (jobid {})".format(jobid))
        self.submitted.append(jobid)
        self.session.deleteJobTemplate(jt)

        thread = threading.Thread(
            target=self._wait_for_job,
            args=(job, jobid, callback, error_callback,
                jobscript))
        thread.daemon = True
        thread.start()
        self.threads.append(thread)

        submit_callback(job)

    def shutdown(self):
        super().shutdown()
        self.session.exit()

    def _wait_for_job(
        self, job, jobid, callback, error_callback,
        jobscript):
        import drmaa
        try:
            retval = self.session.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        except drmaa.errors.InternalException as e:
            print_exception(WorkflowError("DRMAA Error: {}".format(e)), self.workflow.linemaps)
            os.remove(jobscript)
            error_callback(job)
            return
        os.remove(jobscript)
        if retval.exitStatus == 0:
            self.finish_job(job)
            callback(job)
        else:
            self.print_job_error(job)
            print_exception(
                ClusterJobException(job, self.dag.jobid(job), jobscript),
                self.workflow.linemaps)
            error_callback(job)


def run_wrapper(
    run,
    input,
    output,
    params,
    wildcards,
    threads,
    resources,
    log,
    version,
    benchmark,
    benchmark_repeats,
    linemaps
):
    """
    Wrapper around the run method that handles directory creation and
    output file deletion on error.

    Arguments
    run       -- the run method
    input     -- list of input files
    output    -- list of output files
    wildcards -- so far processed wildcards
    threads   -- usable threads
    log       -- path to log file
    """

    if log is None:
        log = Unformattable(errormsg="log used but undefined")
    try:
        runs = 1 if benchmark is None else benchmark_repeats
        wallclock = []
        for i in range(runs):
            w = time.time()
            # execute the actual run method.
            run(input, output, params, wildcards, threads, resources, log, version)
            w = time.time() - w
            wallclock.append(w)

    except (KeyboardInterrupt, SystemExit) as e:
        # re-raise the keyboard interrupt in order to record an error in the scheduler but ignore it
        raise e
    except (Exception, BaseException) as ex:
        # this ensures that exception can be re-raised in the parent thread
        lineno, file = get_exception_origin(ex, linemaps)
        raise RuleException(format_error(
            ex, lineno, linemaps=linemaps, snakefile=file,
            show_traceback=True))

    if benchmark is not None:
        try:
            with open(benchmark, "w") as f:
                json.dump({
                    name: {
                        "s": times,
                        "h:m:s": [str(datetime.timedelta(seconds=t)) for t in times]
                    }
                    for name, times in zip(
                        "wall_clock_times".split(),
                        [wallclock]
                    )
                }, f, indent=4)
        except (Exception, BaseException) as ex:
            raise WorkflowError(ex)
