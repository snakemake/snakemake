__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import shutil
import textwrap
import time
import tarfile
from collections import defaultdict, Counter
from itertools import chain, combinations, filterfalse, product, groupby
from functools import partial, lru_cache
from inspect import isfunction, ismethod
from operator import itemgetter, attrgetter
from pathlib import Path
import subprocess
import uuid

from snakemake.io import IOFile, _IOFile, PeriodicityDetector, wait_for_files, is_flagged, contains_wildcard
from snakemake.jobs import Job, Reason, GroupJob
from snakemake.exceptions import RuleException, MissingInputException
from snakemake.exceptions import MissingRuleException, AmbiguousRuleException
from snakemake.exceptions import CyclicGraphException, MissingOutputException
from snakemake.exceptions import IncompleteFilesException, ImproperOutputException
from snakemake.exceptions import PeriodicWildcardError, WildcardError
from snakemake.exceptions import RemoteFileException, WorkflowError, ChildIOException
from snakemake.exceptions import UnexpectedOutputException, InputFunctionException
from snakemake.logging import logger
from snakemake.common import DYNAMIC_FILL
from snakemake import conda, singularity
from snakemake import utils
from snakemake.output_index import OutputIndex


class DAG:
    """Directed acyclic graph of jobs."""
    def __init__(self,
                 workflow,
                 rules=None,
                 dryrun=False,
                 targetfiles=None,
                 targetrules=None,
                 forceall=False,
                 forcerules=None,
                 forcefiles=None,
                 priorityfiles=None,
                 priorityrules=None,
                 untilfiles=None,
                 untilrules=None,
                 omitfiles=None,
                 omitrules=None,
                 ignore_ambiguity=False,
                 force_incomplete=False,
                 ignore_incomplete=False,
                 notemp=False,
                 keep_remote_local=False):

        self.dryrun = dryrun
        self.dependencies = defaultdict(partial(defaultdict, set))
        self.depending = defaultdict(partial(defaultdict, set))
        self._needrun = set()
        self._priority = dict()
        self._reason = defaultdict(Reason)
        self._finished = set()
        self._dynamic = set()
        self._len = 0
        self.workflow = workflow
        self.rules = set(rules)
        self.ignore_ambiguity = ignore_ambiguity
        self.targetfiles = targetfiles
        self.targetrules = targetrules
        self.priorityfiles = priorityfiles
        self.priorityrules = priorityrules
        self.targetjobs = set()
        self.prioritytargetjobs = set()
        self._ready_jobs = set()
        self.notemp = notemp
        self.keep_remote_local = keep_remote_local
        self._jobid = dict()
        self.job_cache = dict()
        self.conda_envs = dict()
        self.singularity_imgs = dict()
        self._progress = 0
        self._group = dict()

        self.forcerules = set()
        self.forcefiles = set()
        self.untilrules = set()
        self.untilfiles = set()
        self.omitrules = set()
        self.omitfiles = set()
        self.updated_subworkflow_files = set()
        if forceall:
            self.forcerules.update(self.rules)
        elif forcerules:
            self.forcerules.update(forcerules)
        if forcefiles:
            self.forcefiles.update(forcefiles)
        if untilrules:  # keep only the rule names
            self.untilrules.update(set(rule.name for rule in untilrules))
        if untilfiles:
            self.untilfiles.update(untilfiles)
        if omitrules:
            self.omitrules.update(set(rule.name for rule in omitrules))
        if omitfiles:
            self.omitfiles.update(omitfiles)

        self.omitforce = set()

        self.force_incomplete = force_incomplete
        self.ignore_incomplete = ignore_incomplete

        self.periodic_wildcard_detector = PeriodicityDetector()

        self.update_output_index()

    def init(self, progress=False):
        """ Initialise the DAG. """
        for job in map(self.rule2job, self.targetrules):
            job = self.update([job], progress=progress)
            self.targetjobs.add(job)

        for file in self.targetfiles:
            job = self.update(self.file2jobs(file), file=file, progress=progress)
            self.targetjobs.add(job)

        self.cleanup()

        self.update_needrun()
        self.set_until_jobs()
        self.delete_omitfrom_jobs()
        self.update_jobids()

        # Check if there are files/dirs that are children of other outputs.
        allfiles = {}

        for job in self.jobs:
            # This is to account also for targets of symlinks
            allfiles.update({f(x):"input" for x in job.input for f in (os.path.abspath, os.path.realpath)})
            allfiles.update({f(x):"output" for x in job.output for f in (os.path.abspath, os.path.realpath)})

        sortedfiles = sorted(allfiles.keys())
        for i in range(len(sortedfiles)-1):
            if allfiles[sortedfiles[i]] == "output":
                if os.path.commonpath([sortedfiles[i]]) == os.path.commonpath([sortedfiles[i], sortedfiles[i+1]]):
                    raise ChildIOException(parent = sortedfiles[i], child = sortedfiles[i+1])

        # check if remaining jobs are valid
        for i, job in enumerate(self.jobs):
            job.is_valid()

    def update_jobids(self):
        for job in self.jobs:
            if job not in self._jobid:
                self._jobid[job] = len(self._jobid)

    def cleanup(self):
        self.job_cache.clear()
        final_jobs = set(self.jobs)
        todelete = [job for job in self.dependencies if job not in final_jobs]
        for job in todelete:
            del self.dependencies[job]
            try:
                del self.depending[job]
            except KeyError:
                pass

    def create_conda_envs(self,
                          dryrun=False,
                          forceall=False,
                          init_only=False,
                          quiet=False):
        conda.check_conda()
        # First deduplicate based on job.conda_env_file
        jobs = self.jobs if forceall else self.needrun_jobs
        env_set = {(job.conda_env_file, job.singularity_img_url)
                   for job in jobs if job.conda_env_file}
        # Then based on md5sum values
        self.conda_envs = dict()
        for (env_file, simg_url) in env_set:
            simg = None
            if simg_url:
                assert simg_url in self.singularity_imgs, "bug: must first pull singularity images"
                simg = self.singularity_imgs[simg_url]
            env = conda.Env(env_file, self, singularity_img=simg)
            self.conda_envs[(env_file, simg_url)] = env

        if not init_only:
            for env in self.conda_envs.values():
                if not dryrun or not quiet:
                    env.create(dryrun)

    def pull_singularity_imgs(self, dryrun=False, forceall=False, quiet=False):
        # First deduplicate based on job.conda_env_file
        jobs = self.jobs if forceall else self.needrun_jobs
        img_set = {job.singularity_img_url for job in jobs
                   if job.singularity_img_url}

        for img_url in img_set:
            img = singularity.Image(img_url, self)
            if not dryrun or not quiet:
                img.pull(dryrun)
            self.singularity_imgs[img_url] = img


    def update_output_index(self):
        """Update the OutputIndex."""
        self.output_index = OutputIndex(self.rules)

    def check_incomplete(self):
        """Check if any output files are incomplete. This is done by looking up
        markers in the persistence module."""
        if not self.ignore_incomplete:
            incomplete = self.incomplete_files
            if incomplete:
                if self.force_incomplete:
                    logger.debug("Forcing incomplete files:")
                    logger.debug("\t" + "\n\t".join(incomplete))
                    self.forcefiles.update(incomplete)
                else:
                    raise IncompleteFilesException(incomplete)

    def incomplete_external_jobid(self, job):
        """Return the external jobid of the job if it is marked as incomplete.

        Returns None, if job is not incomplete, or if no external jobid has been
        registered or if force_incomplete is True.
        """
        if self.force_incomplete:
            return None
        jobids = self.workflow.persistence.external_jobids(job)
        if len(jobids) == 1:
            return jobids[0]
        elif len(jobids) > 1:
            raise WorkflowError(
                "Multiple different external jobids registered "
                "for output files of incomplete job {} ({}). This job "
                "cannot be resumed. Execute Snakemake with --rerun-incomplete "
                "to fix this issue.".format(job.jobid, jobids))

    def check_dynamic(self):
        """Check dynamic output and update downstream rules if necessary."""
        for job in filter(
                lambda job: (job.dynamic_output and not self.needrun(job)),
                self.jobs):
            self.update_dynamic(job)
        self.postprocess()

    @property
    def dynamic_output_jobs(self):
        """Iterate over all jobs with dynamic output files."""
        return (job for job in self.jobs if job.dynamic_output)

    @property
    def jobs(self):
        """ All jobs in the DAG. """
        for job in self.bfs(self.dependencies, *self.targetjobs):
            yield job

    @property
    def needrun_jobs(self):
        """ Jobs that need to be executed. """
        for job in filter(self.needrun,
                          self.bfs(self.dependencies,
                                   *self.targetjobs,
                                   stop=self.noneedrun_finished)):
            yield job

    @property
    def local_needrun_jobs(self):
        """Iterate over all jobs that need to be run and are marked as local."""
        return filter(lambda job: job.is_local,
                      self.needrun_jobs)

    @property
    def finished_jobs(self):
        """ Iterate over all jobs that have been finished."""
        for job in filter(self.finished, self.bfs(self.dependencies,
                                                  *self.targetjobs)):
            yield job

    @property
    def ready_jobs(self):
        """Jobs that are ready to execute."""
        return self._ready_jobs

    def needrun(self, job):
        """Return whether a given job needs to be executed."""
        return job in self._needrun

    def priority(self, job):
        """Return priority of given job."""
        return self._priority[job]

    def noneedrun_finished(self, job):
        """
        Return whether a given job is finished or was not
        required to run at all.
        """
        return not self.needrun(job) or self.finished(job)

    def reason(self, job):
        """ Return the reason of the job execution. """
        return self._reason[job]

    def finished(self, job):
        """ Return whether a job is finished. """
        return job in self._finished

    def dynamic(self, job):
        """
        Return whether a job is dynamic (i.e. it is only a placeholder
        for those that are created after the job with dynamic output has
        finished.
        """
        if job.is_group():
            for j in job:
                if j in self._dynamic:
                    return True
        else:
            return job in self._dynamic

    def requested_files(self, job):
        """Return the files a job requests."""
        return set(*self.depending[job].values())

    @property
    def incomplete_files(self):
        """Return list of incomplete files."""
        return list(chain(*(job.output
                            for job in
                            filter(self.workflow.persistence.incomplete,
                                   filterfalse(self.needrun, self.jobs)))))

    @property
    def newversion_files(self):
        """Return list of files where the current version is newer than the
        recorded version.
        """
        return list(chain(*(job.output
                            for job in
                            filter(self.workflow.persistence.newversion,
                                   self.jobs))))

    def missing_temp(self, job):
        """
        Return whether a temp file that is input of the given job is missing.
        """
        for job_, files in self.depending[job].items():
            if self.needrun(job_) and any(not f.exists for f in files):
                return True
        return False

    def check_and_touch_output(self,
                               job,
                               wait=3,
                               ignore_missing_output=False,
                               no_touch=False,
                               force_stay_on_remote=False):
        """ Raise exception if output files of job are missing. """
        expanded_output = [job.shadowed_path(path) for path in job.expanded_output]
        if job.benchmark:
            expanded_output.append(job.benchmark)

        if not ignore_missing_output:
            try:
                wait_for_files(expanded_output,
                               latency_wait=wait,
                               force_stay_on_remote=force_stay_on_remote)
            except IOError as e:
                raise MissingOutputException(str(e) + "\nThis might be due to "
                "filesystem latency. If that is the case, consider to increase the "
                "wait time with --latency-wait.", rule=job.rule)

        # Ensure that outputs are of the correct type (those flagged with directory()
        # are directories and not files and vice versa).
        for f in expanded_output:
            if (f.is_directory and not os.path.isdir(f)) or (os.path.isdir(f) and not f.is_directory):
                raise ImproperOutputException(job.rule, [f])

        #It is possible, due to archive expansion or cluster clock skew, that
        #the files appear older than the input.  But we know they must be new,
        #so touch them to update timestamps. This also serves to touch outputs
        #when using the --touch flag.
        #Note that if the input files somehow have a future date then this will
        #not currently be spotted and the job will always be re-run.
        if not no_touch:
            for f in expanded_output:
                # This won't create normal files if missing, but will create
                # the flag file for directories.
                if f.exists_local:
                    f.touch()

    def unshadow_output(self, job, only_log=False):
        """ Move files from shadow directory to real output paths. """
        if not job.shadow_dir or not job.expanded_output:
            return

        files = job.log if only_log else chain(job.expanded_output, job.log)

        for real_output in files:
            shadow_output = job.shadowed_path(real_output).file
            # Remake absolute symlinks as relative
            if os.path.islink(shadow_output):
                dest = os.readlink(shadow_output)
                if os.path.isabs(dest):
                    rel_dest = os.path.relpath(dest, job.shadow_dir)
                    os.remove(shadow_output)
                    os.symlink(rel_dest, shadow_output)

            if os.path.realpath(shadow_output) == os.path.realpath(
                    real_output):
                continue
            logger.debug("Moving shadow output {} to destination {}".format(
                shadow_output, real_output))
            shutil.move(shadow_output, real_output)
        shutil.rmtree(job.shadow_dir)

    def check_periodic_wildcards(self, job):
        """ Raise an exception if a wildcard of the given job appears to be periodic,
        indicating a cyclic dependency. """
        for wildcard, value in job.wildcards_dict.items():
            periodic_substring = self.periodic_wildcard_detector.is_periodic(
                value)
            if periodic_substring is not None:
                raise PeriodicWildcardError(
                    "The value {} in wildcard {} is periodically repeated ({}). "
                    "This would lead to an infinite recursion. "
                    "To avoid this, e.g. restrict the wildcards in this rule to certain values.".format(
                        periodic_substring, wildcard, value),
                    rule=job.rule)

    def handle_protected(self, job):
        """ Write-protect output files that are marked with protected(). """
        for f in job.expanded_output:
            if f in job.protected_output:
                logger.info("Write-protecting output file {}.".format(f))
                f.protect()

    def handle_touch(self, job):
        """ Touches those output files that are marked for touching. """
        for f in job.expanded_output:
            if f in job.touch_output:
                f = job.shadowed_path(f)
                logger.info("Touching output file {}.".format(f))
                f.touch_or_create()
                assert os.path.exists(f)

    def temp_input(self, job):
        for job_, files in self.dependencies[job].items():
            for f in filter(job_.temp_output.__contains__, files):
                yield f

    def temp_size(self, job):
        """Return the total size of temporary input files of the job.
        If none, return 0.
        """
        return sum(f.size for f in self.temp_input(job))

    def handle_temp(self, job):
        """ Remove temp files if they are no longer needed. Update temp_mtimes. """
        if self.notemp:
            return

        is_temp = lambda f: is_flagged(f, "temp")

        # handle temp input
        needed = lambda job_, f: any(
            f in files for j, files in self.depending[job_].items()
            if not self.finished(j) and self.needrun(j) and j != job)

        def unneeded_files():
            # temp input
            for job_, files in self.dependencies[job].items():
                tempfiles = set(f for f in job_.expanded_output if is_temp(f))
                yield from filterfalse(partial(needed, job_), tempfiles & files)

            # temp output
            if not job.dynamic_output:
                tempfiles = (f for f in job.expanded_output
                               if is_temp(f) and f not in self.targetfiles)
                yield from filterfalse(partial(needed, job), tempfiles)

        for f in unneeded_files():
            logger.info("Removing temporary output file {}.".format(f))
            f.remove(remove_non_empty_dir=True)

    def handle_log(self, job, upload_remote=True):
        for f in job.log:
            if not f.exists_local:
                # If log file was not created during job, create an empty one.
                f.touch_or_create()
            if upload_remote and f.is_remote and not f.should_stay_on_remote:
                f.upload_to_remote()
                if not f.exists_remote:
                    raise RemoteFileException(
                        "The file upload was attempted, but it does not "
                        "exist on remote. Check that your credentials have "
                        "read AND write permissions.")

    def handle_remote(self, job, upload=True):
        """ Remove local files if they are no longer needed and upload. """
        if upload:
            # handle output files
            files = list(job.expanded_output)
            if job.benchmark:
                files.append(job.benchmark)
            for f in files:
                if f.is_remote and not f.should_stay_on_remote:
                    f.upload_to_remote()
                    remote_mtime = f.mtime
                    # immediately force local mtime to match remote,
                    # since conversions from S3 headers are not 100% reliable
                    # without this, newness comparisons may fail down the line
                    f.touch(times=(remote_mtime, remote_mtime))

                    if not f.exists_remote:
                        raise RemoteFileException(
                            "The file upload was attempted, but it does not "
                            "exist on remote. Check that your credentials have "
                            "read AND write permissions.")

        if not self.keep_remote_local:
            # handle input files
            needed = lambda job_, f: any(
                f in files for j, files in self.depending[job_].items()
                if not self.finished(j) and self.needrun(j) and j != job)

            def unneeded_files():
                putative = lambda f: f.is_remote and not f.protected and not f.should_keep_local
                generated_input = set()
                for job_, files in self.dependencies[job].items():
                    generated_input |= files
                    for f in filter(putative, files):
                        if not needed(job_, f):
                            yield f
                for f, f_ in zip(job.output, job.rule.output):
                    if putative(f) and not needed(job, f) and not f in self.targetfiles:
                        if f in job.dynamic_output:
                            for f_ in job.expand_dynamic(f_):
                                yield f_
                        else:
                            yield f
                for f in filter(putative, job.input):
                    # TODO what about remote inputs that are used by multiple jobs?
                    if f not in generated_input:
                        yield f

            for f in unneeded_files():
                if f.exists_local:
                    logger.info("Removing local output file: {}".format(f))
                    f.remove()

            job.rmdir_empty_remote_dirs()

    def jobid(self, job):
        """Return job id of given job."""
        return self._jobid[job]

    def update(self, jobs, file=None, visited=None, skip_until_dynamic=False, progress=False):
        """ Update the DAG by adding given jobs and their dependencies. """
        if visited is None:
            visited = set()
        producer = None
        exceptions = list()
        jobs = sorted(jobs, reverse=not self.ignore_ambiguity)
        cycles = list()


        for job in jobs:
            logger.dag_debug(dict(status="candidate", job=job))
            if file in job.input:
                cycles.append(job)
                continue
            if job in visited:
                cycles.append(job)
                continue
            try:
                self.check_periodic_wildcards(job)
                self.update_(job,
                             visited=set(visited),
                             skip_until_dynamic=skip_until_dynamic,
                             progress=progress)
                # TODO this might fail if a rule discarded here is needed
                # elsewhere
                if producer:
                    if job < producer or self.ignore_ambiguity:
                        break
                    elif producer is not None:
                        raise AmbiguousRuleException(file, job, producer)
                producer = job
            except (MissingInputException, CyclicGraphException,
                    PeriodicWildcardError) as ex:
                exceptions.append(ex)
            except RecursionError as e:
                raise WorkflowError(
                    e, "If building the DAG exceeds the recursion limit, "
                    "this is likely due to a cyclic dependency."
                    "E.g. you might have a sequence of rules that "
                    "can generate their own input. Try to make "
                    "the output files more specific. "
                    "A common pattern is to have different prefixes "
                    "in the output files of different rules." +
                    "\nProblematic file pattern: {}".format(file) if file else
                    "")
        if producer is None:
            if cycles:
                job = cycles[0]
                raise CyclicGraphException(job.rule, file, rule=job.rule)
            if exceptions:
                raise exceptions[0]
        else:
            logger.dag_debug(dict(status="selected", job=producer))

        n = len(self.dependencies)
        if progress and n % 1000 == 0 and n and self._progress != n:
            logger.info("Processed {} potential jobs.".format(n))
            self._progress = n

        return producer

    def update_(self, job, visited=None, skip_until_dynamic=False, progress=False):
        """ Update the DAG by adding the given job and its dependencies. """
        if job in self.dependencies:
            return
        if visited is None:
            visited = set()
        visited.add(job)
        dependencies = self.dependencies[job]
        potential_dependencies = self.collect_potential_dependencies(
            job).items()

        skip_until_dynamic = skip_until_dynamic and not job.dynamic_output

        missing_input = job.missing_input
        producer = dict()
        exceptions = dict()
        for file, jobs in potential_dependencies:
            try:
                selected_job = self.update(
                    jobs,
                    file=file,
                    visited=visited,
                    skip_until_dynamic=skip_until_dynamic or file in
                    job.dynamic_input,
                    progress=progress)
                producer[file] = selected_job
            except (MissingInputException, CyclicGraphException,
                    PeriodicWildcardError) as ex:
                if file in missing_input:
                    self.delete_job(job,
                                    recursive=False)  # delete job from tree
                    raise ex

        for file, job_ in producer.items():
            dependencies[job_].add(file)
            self.depending[job_][job].add(file)

        missing_input -= producer.keys()
        if missing_input:
            self.delete_job(job, recursive=False)  # delete job from tree
            raise MissingInputException(job.rule, missing_input)

        if skip_until_dynamic:
            self._dynamic.add(job)

    def update_needrun(self):
        """ Update the information whether a job needs to be executed. """

        def output_mintime(job):
            for job_ in self.bfs(self.depending, job):
                t = job_.output_mintime
                if t:
                    return t

        def needrun(job):
            reason = self.reason(job)
            noinitreason = not reason
            updated_subworkflow_input = self.updated_subworkflow_files.intersection(
                job.input)
            if (job not in self.omitforce and job.rule in self.forcerules or
                    not self.forcefiles.isdisjoint(job.output)):
                reason.forced = True
            elif updated_subworkflow_input:
                reason.updated_input.update(updated_subworkflow_input)
            elif job in self.targetjobs:
                # TODO find a way to handle added/removed input files here?
                if not job.output and not job.benchmark:
                    if job.input:
                        if job.rule.norun:
                            reason.updated_input_run.update(
                                [f for f in job.input if not f.exists])
                        else:
                            reason.nooutput = True
                    else:
                        reason.noio = True
                else:
                    if job.rule in self.targetrules:
                        missing_output = job.missing_output()
                    else:
                        missing_output = job.missing_output(
                            requested=set(chain(*self.depending[job].values(
                            ))) | self.targetfiles)
                    reason.missing_output.update(missing_output)
            if not reason:
                output_mintime_ = output_mintime(job)
                if output_mintime_:
                    updated_input = [
                        f
                        for f in job.input
                        if f.exists and f.is_newer(output_mintime_)
                    ]
                    reason.updated_input.update(updated_input)
            if noinitreason and reason:
                reason.derived = False
            return job

        reason = self.reason
        _needrun = self._needrun
        dependencies = self.dependencies
        depending = self.depending

        _needrun.clear()
        candidates = set(self.jobs)

        queue = list(filter(reason, map(needrun, candidates)))
        visited = set(queue)
        while queue:
            job = queue.pop(0)
            _needrun.add(job)

            for job_, files in dependencies[job].items():
                missing_output = job_.missing_output(requested=files)
                reason(job_).missing_output.update(missing_output)
                if missing_output and not job_ in visited:
                    visited.add(job_)
                    queue.append(job_)

            for job_, files in depending[job].items():
                if job_ in candidates:
                    reason(job_).updated_input_run.update(files)
                    if not job_ in visited:
                        visited.add(job_)
                        queue.append(job_)

        self._len = len(_needrun)

    def in_until(self, job):
        """Return whether given job has been specified via --until."""
        return (job.rule.name in self.untilrules or
                not self.untilfiles.isdisjoint(job.output))

    def in_omitfrom(self, job):
        """Return whether given job has been specified via --omit-from."""
        return (job.rule.name in self.omitrules or
                not self.omitfiles.isdisjoint(job.output))

    def until_jobs(self):
        """Returns a generator of jobs specified by untiljobs."""
        return (job for job in self.jobs if self.in_until(job))

    def omitfrom_jobs(self):
        """Returns a generator of jobs specified by omitfromjobs."""
        return (job for job in self.jobs if self.in_omitfrom(job))

    def downstream_of_omitfrom(self):
        """Returns the downstream of --omit-from rules or files."""
        return filter(lambda job: not self.in_omitfrom(job),
                      self.bfs(self.depending, *self.omitfrom_jobs()))

    def delete_omitfrom_jobs(self):
        """Removes jobs downstream of jobs specified by --omit-from."""
        if not self.omitrules and not self.omitfiles:
            return
        downstream_jobs = list(self.downstream_of_omitfrom()
                               )  # need to cast as list before deleting jobs
        for job in downstream_jobs:
            self.delete_job(job, recursive=False, add_dependencies=True)

    def set_until_jobs(self):
        """Removes jobs downstream of jobs specified by --omit-from."""
        if not self.untilrules and not self.untilfiles:
            return
        self.targetjobs = set(self.until_jobs())

    def update_priority(self):
        """ Update job priorities. """
        prioritized = (
            lambda job: job.rule in self.priorityrules or not self.priorityfiles.isdisjoint(job.output)
        )
        for job in self.needrun_jobs:
            self._priority[job] = job.rule.priority
        for job in self.bfs(self.dependencies,
                            *filter(prioritized, self.needrun_jobs),
                            stop=self.noneedrun_finished):
            self._priority[job] = Job.HIGHEST_PRIORITY

    def update_groups(self):
        groups = dict()
        for job in self.needrun_jobs:
            if job.group is None:
                continue
            stop = lambda j: j.group != job.group
            # BFS into depending jobs if in same group
            # Note: never go up here (into depending), because it may contain
            # jobs that have been sorted out due to e.g. ruleorder.
            group = GroupJob(job.group,
                             self.bfs(self.dependencies, job, stop=stop))

            # merge with previously determined groups if present
            for j in group:
                if j in groups:
                    other = groups[j]
                    other.merge(group)
                    group = other
            # update assignment
            for j in group:
                if j not in groups:
                    groups[j] = group
        self._group = groups

    def update_ready(self, jobs=None):
        """ Update information whether a job is ready to execute.

        Given jobs must be needrun jobs!
        """

        if jobs is None:
            jobs = self.needrun_jobs

        candidate_groups = set()
        for job in jobs:
            if not self.finished(job) and self._ready(job):
                if job.group is None:
                    self._ready_jobs.add(job)
                else:
                    candidate_groups.add(self._group[job])

        self._ready_jobs.update(
            group for group in candidate_groups
            if all(self._ready(job) for job in group))

    def get_jobs_or_groups(self):
        visited_groups = set()
        for job in self.jobs:
            if job.group is None:
                yield job
            else:
                group = self._group[job]
                if group in visited_groups:
                    continue
                visited_groups.add(group)
                yield group

    def close_remote_objects(self):
        """Close all remote objects."""
        for job in self.jobs:
            if not self.needrun(job):
                job.close_remote()

    def postprocess(self):
        """Postprocess the DAG. This has to be invoked after any change to the
        DAG topology."""
        self.update_jobids()
        self.update_needrun()
        self.update_priority()
        self.handle_pipes()
        self.update_groups()
        self.update_ready()
        self.close_remote_objects()

    def handle_pipes(self):
        """Use pipes to determine job groups. Check if every pipe has exactly
           one consumer"""
        for job in self.needrun_jobs:
            candidate_groups = set()
            if job.group is not None:
                candidate_groups.add(job.group)
            all_depending = set()
            has_pipe = False
            for f in job.output:
                if is_flagged(f, "pipe"):
                    if job.is_run:
                        raise WorkflowError("Rule defines pipe output but "
                                            "uses a 'run' directive. This is "
                                            "not possible for technical "
                                            "reasons. Consider using 'shell' or "
                                            "'script'.", rule=job.rule)

                    has_pipe = True
                    depending = [j for j, files in self.depending[job].items()
                                   if f in files]
                    if len(depending) > 1:
                        raise WorkflowError("Output file {} is marked as pipe "
                                            "but more than one job depends on "
                                            "it. Make sure that any pipe "
                                            "output is only consumed by one "
                                            "job".format(f),
                                            rule=job.rule)
                    elif len(depending) == 0:
                        raise WorkflowError("Output file {} is marked as pipe "
                                            "but it has no consumer. This is "
                                            "invalid because it can lead to "
                                            "a dead lock.".format(f),
                                            rule=job.rule)

                    depending = depending[0]

                    if depending.is_run:
                        raise WorkflowError("Rule consumes pipe input but "
                                            "uses a 'run' directive. This is "
                                            "not possible for technical "
                                            "reasons. Consider using 'shell' or "
                                            "'script'.", rule=job.rule)

                    all_depending.add(depending)
                    if depending.group is not None:
                        candidate_groups.add(depending.group)
            if not has_pipe:
                continue

            if len(candidate_groups) > 1:
                raise WorkflowError("An output file is marked as "
                                    "pipe, but consuming jobs "
                                    "are part of conflicting "
                                    "groups.",
                                    rule=job.rule)
            elif candidate_groups:
                # extend the candidate group to all involved jobs
                group = candidate_groups.pop()
            else:
                # generate a random unique group name
                group = uuid.uuid4()
            job.group = group
            for j in all_depending:
                j.group = group

    def _ready(self, job):
        """Return whether the given job is ready to execute."""
        group = self._group.get(job, None)
        if group is None:
            is_external_needrun_dep = self.needrun
        else:
            def is_external_needrun_dep(j):
                g = self._group.get(j, None)
                return self.needrun(j) and (g is None or g != group)

        return self._finished.issuperset(filter(is_external_needrun_dep,
                                                self.dependencies[job]))

    def finish(self, job, update_dynamic=True):
        """Finish a given job (e.g. remove from ready jobs, mark depending jobs
        as ready)."""
        try:
            self._ready_jobs.remove(job)
        except KeyError:
            pass

        if job.is_group():
            jobs = job
        else:
            jobs = [job]

        self._finished.update(jobs)

        # mark depending jobs as ready
        # skip jobs that are marked as until jobs
        self.update_ready(j for job in jobs for j in self.depending[job]
                            if not self.in_until(job) and self.needrun(j))

        for job in jobs:
            if update_dynamic and job.dynamic_output:
                logger.info("Dynamically updating jobs")
                newjob = self.update_dynamic(job)
                if newjob:
                    # simulate that this job ran and was finished before
                    self.omitforce.add(newjob)
                    self._needrun.add(newjob)
                    self._finished.add(newjob)

                    self.postprocess()
                    self.handle_protected(newjob)
                    self.handle_touch(newjob)
                    # add finished jobs to len as they are not counted after new postprocess
                    self._len += len(self._finished)

    def new_job(self, rule, targetfile=None, format_wildcards=None):
        """Create new job for given rule and (optional) targetfile.
        This will reuse existing jobs with the same wildcards."""
        key = (rule, targetfile)
        if key in self.job_cache:
            assert targetfile is not None
            return self.job_cache[key]
        wildcards_dict = rule.get_wildcards(targetfile)
        job = Job(rule, self, wildcards_dict=wildcards_dict,
                  format_wildcards=format_wildcards, targetfile=targetfile)
        for f in job.products:
            self.job_cache[(rule, f)] = job
        return job

    def update_dynamic(self, job):
        """Update the DAG by evaluating the output of the given job that
        contains dynamic output files."""
        dynamic_wildcards = job.dynamic_wildcards
        if not dynamic_wildcards:
            # this happens e.g. in dryrun if output is not yet present
            return

        depending = list(filter(lambda job_: not self.finished(job_), self.bfs(
            self.depending, job)))
        newrule, non_dynamic_wildcards = job.rule.dynamic_branch(
            dynamic_wildcards,
            input=False)
        self.specialize_rule(job.rule, newrule)

        # no targetfile needed for job
        newjob = self.new_job(newrule, format_wildcards=non_dynamic_wildcards)
        self.replace_job(job, newjob)
        for job_ in depending:
            needs_update = any(
                f.get_wildcard_names() & dynamic_wildcards.keys()
                for f in job_.rule.dynamic_input)

            if needs_update:
                newrule_ = job_.rule.dynamic_branch(dynamic_wildcards)
                if newrule_ is not None:
                    self.specialize_rule(job_.rule, newrule_)
                    if not self.dynamic(job_):
                        logger.debug("Updating job {}.".format(job_))
                        newjob_ = self.new_job(
                            newrule_,
                            targetfile=job_.output[0] if job_.output else None)

                        unexpected_output = self.reason(
                            job_).missing_output.intersection(
                                newjob.existing_output)
                        if unexpected_output:
                            logger.warning(
                                "Warning: the following output files of rule {} were not "
                                "present when the DAG was created:\n{}".format(
                                    newjob_.rule, unexpected_output))

                        self.replace_job(job_, newjob_)
        return newjob

    def delete_job(self, job, recursive=True, add_dependencies=False):
        """Delete given job from DAG."""
        if job in self.targetjobs:
            self.targetjobs.remove(job)
        if add_dependencies:
            for _job in self.dependencies[job]:
                self.targetjobs.add(_job)
        for job_ in self.depending[job]:
            del self.dependencies[job_][job]
        del self.depending[job]
        for job_ in self.dependencies[job]:
            depending = self.depending[job_]
            del depending[job]
            if not depending and recursive:
                self.delete_job(job_)
        del self.dependencies[job]
        if job in self._needrun:
            self._len -= 1
            self._needrun.remove(job)
            del self._reason[job]
        if job in self._finished:
            self._finished.remove(job)
        if job in self._dynamic:
            self._dynamic.remove(job)
        if job in self._ready_jobs:
            self._ready_jobs.remove(job)

    def replace_job(self, job, newjob):
        """Replace given job with new job."""
        if job in self.targetjobs:
            self.targetjobs.remove(job)
            self.targetjobs.add(newjob)
        depending = list(self.depending[job].items())
        if self.finished(job):
            self._finished.add(newjob)

        self.delete_job(job)
        self.update([newjob])

        logger.debug("Replace {} with dynamic branch {}".format(job, newjob))
        for job_, files in depending:
            #if not job_.dynamic_input:
            logger.debug("updating depending job {}".format(job_))
            self.dependencies[job_][newjob].update(files)
            self.depending[newjob][job_].update(files)

    def specialize_rule(self, rule, newrule):
        """Specialize the given rule by inserting newrule into the DAG."""
        assert newrule is not None
        self.rules.add(newrule)
        self.update_output_index()

    def collect_potential_dependencies(self, job):
        """Collect all potential dependencies of a job. These might contain
        ambiguities."""
        dependencies = defaultdict(list)
        # use a set to circumvent multiple jobs for the same file
        # if user specified it twice
        file2jobs = self.file2jobs
        for file in job.unique_input:
            # omit the file if it comes from a subworkflow
            if file in job.subworkflow_input:
                continue
            try:
                if file in job.dependencies:
                    jobs = [self.new_job(job.dependencies[file], targetfile=file)]
                else:
                    jobs = file2jobs(file)
                dependencies[file].extend(jobs)
            except MissingRuleException as ex:
                pass
        return dependencies

    def bfs(self, direction, *jobs, stop=lambda job: False):
        """Perform a breadth-first traversal of the DAG."""
        queue = list(jobs)
        visited = set(queue)
        while queue:
            job = queue.pop(0)
            if stop(job):
                # stop criterion reached for this node
                continue
            yield job
            for job_, _ in direction[job].items():
                if not job_ in visited:
                    queue.append(job_)
                    visited.add(job_)

    def level_bfs(self, direction, *jobs, stop=lambda job: False):
        """Perform a breadth-first traversal of the DAG, but also yield the
        level together with each job."""
        queue = [(job, 0) for job in jobs]
        visited = set(jobs)
        while queue:
            job, level = queue.pop(0)
            if stop(job):
                # stop criterion reached for this node
                continue
            yield level, job
            level += 1
            for job_, _ in direction[job].items():
                if not job_ in visited:
                    queue.append((job_, level))
                    visited.add(job_)

    def dfs(self, direction, *jobs, stop=lambda job: False, post=True):
        """Perform depth-first traversal of the DAG."""
        visited = set()

        def _dfs(job):
            """Inner function for DFS traversal."""
            if stop(job):
                return
            if not post:
                yield job
            for job_ in direction[job]:
                if not job_ in visited:
                    visited.add(job_)
                    for j in _dfs(job_):
                        yield j
            if post:
                yield job

        for job in jobs:
            for job_ in self._dfs(direction,
                                  job,
                                  visited,
                                  stop=stop,
                                  post=post):
                yield job_

    def new_wildcards(self, job):
        """Return wildcards that are newly introduced in this job,
        compared to its ancestors."""
        new_wildcards = set(job.wildcards.items())
        for job_ in self.dependencies[job]:
            if not new_wildcards:
                return set()
            for wildcard in job_.wildcards.items():
                new_wildcards.discard(wildcard)
        return new_wildcards

    def rule2job(self, targetrule):
        """Generate a new job from a given rule."""
        if targetrule.has_wildcards():
            raise WorkflowError("Target rules may not contain wildcards. Please specify concrete files or a rule without wildcards.")
        return self.new_job(targetrule)

    def file2jobs(self, targetfile):
        rules = self.output_index.match(targetfile)
        jobs = []
        exceptions = list()
        for rule in rules:
            if rule.is_producer(targetfile):
                try:
                    jobs.append(self.new_job(rule, targetfile=targetfile))
                except InputFunctionException as e:
                    exceptions.append(e)
        if not jobs:
            if exceptions:
                raise exceptions[0]
            raise MissingRuleException(targetfile)
        return jobs

    def rule_dot2(self):
        dag = defaultdict(list)
        visited = set()
        preselect = set()

        def preselect_parents(job):
            for parent in self.depending[job]:
                if parent in preselect:
                    continue
                preselect.add(parent)
                preselect_parents(parent)

        def build_ruledag(job, key=lambda job: job.rule.name):
            if job in visited:
                return
            visited.add(job)
            deps = sorted(self.dependencies[job], key=key)
            deps = [(group[0] if preselect.isdisjoint(group) else
                     preselect.intersection(group).pop())
                    for group in (list(g) for _, g in groupby(deps, key))]
            dag[job].extend(deps)
            preselect_parents(job)
            for dep in deps:
                build_ruledag(dep)

        for job in self.targetjobs:
            build_ruledag(job)

        return self._dot(dag.keys(),
                         print_wildcards=False,
                         print_types=False,
                         dag=dag)

    def rule_dot(self):
        graph = defaultdict(set)
        for job in self.jobs:
            graph[job.rule].update(dep.rule for dep in self.dependencies[job])
        return self._dot(graph)

    def dot(self):
        def node2style(job):
            if not self.needrun(job):
                return "rounded,dashed"
            if self.dynamic(job) or job.dynamic_input:
                return "rounded,dotted"
            return "rounded"

        def format_wildcard(wildcard):
            name, value = wildcard
            if DYNAMIC_FILL in value:
                value = "..."
            return "{}: {}".format(name, value)

        node2rule = lambda job: job.rule
        node2label = lambda job: "\\n".join(chain([
            job.rule.name
        ], sorted(map(format_wildcard, self.new_wildcards(job)))))

        dag = {job: self.dependencies[job] for job in self.jobs}

        return self._dot(dag,
                         node2rule=node2rule,
                         node2style=node2style,
                         node2label=node2label)

    def _dot(self,
             graph,
             node2rule=lambda node: node,
             node2style=lambda node: "rounded",
             node2label=lambda node: node):

        # color rules
        huefactor = 2 / (3 * len(self.rules))
        rulecolor = {
            rule: "{:.2f} 0.6 0.85".format(i * huefactor)
            for i, rule in enumerate(self.rules)
        }

        # markup
        node_markup = '\t{}[label = "{}", color = "{}", style="{}"];'.format
        edge_markup = "\t{} -> {}".format

        # node ids
        ids = {node: i for i, node in enumerate(graph)}

        # calculate nodes
        nodes = [node_markup(ids[node], node2label(node),
                             rulecolor[node2rule(node)], node2style(node))
                 for node in graph]
        # calculate edges
        edges = [edge_markup(ids[dep], ids[node])
                 for node, deps in graph.items() for dep in deps]

        return textwrap.dedent("""\
            digraph snakemake_dag {{
                graph[bgcolor=white, margin=0];
                node[shape=box, style=rounded, fontname=sans, \
                fontsize=10, penwidth=2];
                edge[penwidth=2, color=grey];
            {items}
            }}\
            """).format(items="\n".join(nodes + edges))

    def summary(self, detailed=False):
        if detailed:
            yield "output_file\tdate\trule\tversion\tlog-file(s)\tinput-file(s)\tshellcmd\tstatus\tplan"
        else:
            yield "output_file\tdate\trule\tversion\tlog-file(s)\tstatus\tplan"

        for job in self.jobs:
            output = job.rule.output if self.dynamic(
                job) else job.expanded_output
            for f in output:
                rule = self.workflow.persistence.rule(f)
                rule = "-" if rule is None else rule

                version = self.workflow.persistence.version(f)
                version = "-" if version is None else str(version)

                date = time.ctime(f.mtime) if f.exists else "-"

                pending = "update pending" if self.reason(job) else "no update"

                log = self.workflow.persistence.log(f)
                log = "-" if log is None else ",".join(log)

                input = self.workflow.persistence.input(f)
                input = "-" if input is None else ",".join(input)

                shellcmd = self.workflow.persistence.shellcmd(f)
                shellcmd = "-" if shellcmd is None else shellcmd
                # remove new line characters, leading and trailing whitespace
                shellcmd = shellcmd.strip().replace("\n", "; ")

                status = "ok"
                if not f.exists:
                    status = "missing"
                elif self.reason(job).updated_input:
                    status = "updated input files"
                elif self.workflow.persistence.version_changed(job, file=f):
                    status = "version changed to {}".format(job.rule.version)
                elif self.workflow.persistence.code_changed(job, file=f):
                    status = "rule implementation changed"
                elif self.workflow.persistence.input_changed(job, file=f):
                    status = "set of input files changed"
                elif self.workflow.persistence.params_changed(job, file=f):
                    status = "params changed"
                if detailed:
                    yield "\t".join((f, date, rule, version, log, input, shellcmd,
                                     status, pending))
                else:
                    yield "\t".join((f, date, rule, version, log, status, pending))

    def archive(self, path):
        """Archives workflow such that it can be re-run on a different system.

        Archiving includes git versioned files (i.e. Snakefiles, config files, ...),
        ancestral input files and conda environments.
        """
        if path.endswith(".tar"):
            mode = "x"
        elif path.endswith("tar.bz2"):
            mode = "x:bz2"
        elif path.endswith("tar.xz"):
            mode = "x:xz"
        elif path.endswith("tar.gz"):
            mode = "x:gz"
        else:
            raise WorkflowError("Unsupported archive format "
                                "(supported: .tar, .tar.gz, .tar.bz2, .tar.xz)")
        if os.path.exists(path):
            raise WorkflowError("Archive already exists:\n" + path)

        self.create_conda_envs(forceall=True)

        try:
            workdir = Path(os.path.abspath(os.getcwd()))
            with tarfile.open(path, mode=mode, dereference=True) as archive:
                archived = set()

                def add(path):
                    if workdir not in Path(os.path.abspath(path)).parents:
                        logger.warning("Path {} cannot be archived: "
                                       "not within working directory.".format(path))
                    else:
                        f = os.path.relpath(path)
                        if f not in archived:
                            archive.add(f)
                            archived.add(f)
                            logger.info("archived " + f)

                logger.info("Archiving snakefiles, scripts and files under "
                            "version control...")
                for f in self.workflow.get_sources():
                    add(f)

                logger.info("Archiving external input files...")
                for job in self.jobs:
                    # input files
                    for f in job.input:
                        if not any(f in files for files in self.dependencies[job].values()):
                            # this is an input file that is not created by any job
                            add(f)

                logger.info("Archiving conda environments...")
                envs = set()
                for job in self.jobs:
                    if job.conda_env_file:
                        env_archive = job.archive_conda_env()
                        envs.add(env_archive)
                for env in envs:
                    add(env)

        except (Exception, BaseException) as e:
            os.remove(path)
            raise e


    def clean(self, only_temp=False, dryrun=False):
        """Removes files generated by the workflow.
        """
        if dryrun:
            logger.info("Would delete the following files:")
        else:
            logger.info("Deleting the following files:")

        for job in self.jobs:
            for f in job.output:
                if not only_temp or is_flagged(f, "temp"):
                    # The reason for the second check is that dangling
                    # symlinks fail f.exists.
                    if f.exists or os.path.islink(f):
                        logger.info(f)
                        if not dryrun:
                            # Remove non-empty dirs if flagged as temp()
                            f.remove(remove_non_empty_dir=only_temp)

    def list_untracked(self):
        """List files in the workdir that are not in the dag.
        """
        used_files = set()
        files_in_cwd = set()
        for job in self.jobs:
            used_files.update(os.path.relpath(file) for file in chain(job.local_input, job.local_output, job.log))
        for root, dirs, files in os.walk(os.getcwd()):
            # Ignore hidden files and don't traverse into hidden dirs
            files_in_cwd.update([os.path.relpath(os.path.join(root, f)) for f in files if not f[0] == '.'])
            dirs[:] = [d for d in dirs if not d[0] == '.']
        for f in sorted(list(files_in_cwd - used_files)):
            logger.info(f)

    def d3dag(self, max_jobs=10000):
        def node(job):
            jobid = self.jobid(job)
            return {
                "id": jobid,
                "value": {
                    "jobid": jobid,
                    "label": job.rule.name,
                    "rule": job.rule.name
                }
            }

        def edge(a, b):
            return {"u": self.jobid(a), "v": self.jobid(b)}

        jobs = list(self.jobs)

        if len(jobs) > max_jobs:
            logger.info(
                "Job-DAG is too large for visualization (>{} jobs).".format(
                    max_jobs))
        else:
            logger.d3dag(nodes=[node(job) for job in jobs],
                         edges=[edge(dep, job)
                                for job in jobs for dep in self.dependencies[
                                    job] if self.needrun(dep)])

    def stats(self):
        rules = Counter()
        rules.update(job.rule for job in self.needrun_jobs)
        rules.update(job.rule for job in self.finished_jobs)
        yield "Job counts:"
        yield "\tcount\tjobs"
        for rule, count in sorted(rules.most_common(),
                                  key=lambda item: item[0].name):
            yield "\t{}\t{}".format(count, rule)
        yield "\t{}".format(len(self))

    def __str__(self):
        return self.dot()

    def __len__(self):
        return self._len
