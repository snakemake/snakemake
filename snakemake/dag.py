__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import shutil
import textwrap
import time
from collections import defaultdict, Counter
from itertools import chain, combinations, filterfalse, product, groupby
from functools import partial, lru_cache
from inspect import isfunction, ismethod
from operator import itemgetter, attrgetter

from snakemake.io import IOFile, _IOFile, PeriodicityDetector, wait_for_files, is_flagged, contains_wildcard
from snakemake.jobs import Job, Reason
from snakemake.exceptions import RuleException, MissingInputException
from snakemake.exceptions import MissingRuleException, AmbiguousRuleException
from snakemake.exceptions import CyclicGraphException, MissingOutputException
from snakemake.exceptions import IncompleteFilesException
from snakemake.exceptions import PeriodicWildcardError
from snakemake.exceptions import RemoteFileException, WorkflowError
from snakemake.exceptions import UnexpectedOutputException, InputFunctionException
from snakemake.logging import logger
from snakemake.output_index import OutputIndex
from snakemake.common import DYNAMIC_FILL

# Workaround for Py <3.5 prior to existence of RecursionError
try:
    RecursionError
except NameError:
    RecursionError = RuntimeError


class DAG:
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
                 notemp=False):

        self.dryrun = dryrun
        self.dependencies = defaultdict(partial(defaultdict, set))
        self.depending = defaultdict(partial(defaultdict, set))
        self._needrun = set()
        self._priority = dict()
        self._downstream_size = dict()
        self._temp_input_count = dict()
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
        self._jobid = dict()

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

    def init(self):
        """ Initialise the DAG. """
        for job in map(self.rule2job, self.targetrules):
            job = self.update([job])
            self.targetjobs.add(job)

        for file in self.targetfiles:
            job = self.update(self.file2jobs(file), file=file)
            self.targetjobs.add(job)

        self.update_needrun()
        self.set_until_jobs()
        self.delete_omitfrom_jobs()

    def update_output_index(self):
        self.output_index = OutputIndex(self.rules)

    def check_incomplete(self):
        if not self.ignore_incomplete:
            incomplete = self.incomplete_files
            if incomplete:
                if self.force_incomplete:
                    logger.debug("Forcing incomplete files:")
                    logger.debug("\t" + "\n\t".join(incomplete))
                    self.forcefiles.update(incomplete)
                else:
                    raise IncompleteFilesException(incomplete)

    def check_dynamic(self):
        for job in filter(
                lambda job: (job.dynamic_output and not self.needrun(job)),
                self.jobs):
            self.update_dynamic(job)

    @property
    def dynamic_output_jobs(self):
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
        return filter(lambda job: self.workflow.is_local(job.rule),
                      self.needrun_jobs)

    @property
    def finished_jobs(self):
        """ Jobs that have been executed. """
        for job in filter(self.finished, self.bfs(self.dependencies,
                                                  *self.targetjobs)):
            yield job

    @property
    def ready_jobs(self):
        """ Jobs that are ready to execute. """
        return self._ready_jobs

    def ready(self, job):
        """ Return whether a given job is ready to execute. """
        return job in self._ready_jobs

    def needrun(self, job):
        """ Return whether a given job needs to be executed. """
        return job in self._needrun

    def priority(self, job):
        return self._priority[job]

    def downstream_size(self, job):
        return self._downstream_size[job]

    def temp_input_count(self, job):
        return self._temp_input_count[job]

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
        return job in self._dynamic

    def requested_files(self, job):
        """ Return the files a job requests. """
        return set(*self.depending[job].values())

    @property
    def incomplete_files(self):
        return list(chain(*(job.output
                            for job in
                            filter(self.workflow.persistence.incomplete,
                                   filterfalse(self.needrun, self.jobs)))))

    @property
    def newversion_files(self):
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

    def check_output(self, job, wait=3):
        """ Raise exception if output files of job are missing. """
        try:
            wait_for_files(job.expanded_shadowed_output, latency_wait=wait)
        except IOError as e:
            raise MissingOutputException(str(e), rule=job.rule)

        input_maxtime = job.input_maxtime
        if input_maxtime is not None:
            output_mintime = job.output_mintime_local
            if output_mintime is not None and output_mintime < input_maxtime:
                raise RuleException(
                    "Output files {} are older than input "
                    "files. Did you extract an archive? Make sure that output "
                    "files have a more recent modification date than the "
                    "archive, e.g. by using 'touch'.".format(", ".join(
                        job.expanded_output)),
                    rule=job.rule)

    def unshadow_output(self, job):
        """ Move files from shadow directory to real output paths. """
        if not job.shadow_dir or not job.expanded_output:
            return
        cwd = os.getcwd()
        for real_output in job.expanded_output:
            shadow_output = os.path.join(job.shadow_dir, real_output)
            if os.path.realpath(shadow_output) == os.path.realpath(
                    real_output):
                continue
            logger.info("Moving shadow output {} to destination {}".format(
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
                logger.info("Touching output file {}.".format(f))
                f.touch_or_create()

    def temp_input(self, job):
        for job_, files in self.dependencies[job].items():
            for f in filter(job_.temp_output.__contains__, files):
                yield f

    def handle_temp(self, job):
        """ Remove temp files if they are no longer needed. """
        if self.notemp:
            return

        needed = lambda job_, f: any(
            f in files for j, files in self.depending[job_].items()
            if not self.finished(j) and self.needrun(j) and j != job)

        def unneeded_files():
            for job_, files in self.dependencies[job].items():
                for f in job_.temp_output & files:
                    if not needed(job_, f):
                        yield f
            for f in filterfalse(partial(needed, job), job.temp_output):
                if not f in self.targetfiles:
                    yield f

        for f in unneeded_files():
            logger.info("Removing temporary output file {}.".format(f))
            f.remove(remove_non_empty_dir=True)

    def handle_remote(self, job):
        """ Remove local files if they are no longer needed, and upload to S3. """

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
            for f in filter(putative, job.output):
                if not needed(job, f) and not f in self.targetfiles:
                    for f_ in job.expand_dynamic(f):
                        yield f
            for f in filter(putative, job.input):
                # TODO what about remote inputs that are used by multiple jobs?
                if f not in generated_input:
                    yield f

        for f in job.expanded_output:
            if f.is_remote:
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

        for f in unneeded_files():
            logger.info("Removing local output file: {}".format(f))
            f.remove()

        job.rmdir_empty_remote_dirs()

    def jobid(self, job):
        if job not in self._jobid:
            self._jobid[job] = len(self._jobid)
        return self._jobid[job]

    def update(self, jobs, file=None, visited=None, skip_until_dynamic=False):
        """ Update the DAG by adding given jobs and their dependencies. """
        if visited is None:
            visited = set()
        producer = None
        exceptions = list()
        jobs = sorted(jobs, reverse=not self.ignore_ambiguity)
        cycles = list()

        for job in jobs:
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
                             skip_until_dynamic=skip_until_dynamic)
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
        return producer

    def update_(self, job, visited=None, skip_until_dynamic=False):
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
                producer[file] = self.update(
                    jobs,
                    file=file,
                    visited=visited,
                    skip_until_dynamic=skip_until_dynamic or file in
                    job.dynamic_input)
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
        return (job.rule.name in self.untilrules or
                not self.untilfiles.isdisjoint(job.output))

    def in_omitfrom(self, job):
        return (job.rule.name in self.omitrules or
                not self.omitfiles.isdisjoint(job.output))

    def until_jobs(self):
        'Returns a generator of jobs specified by untiljobs'
        return (job for job in self.jobs if self.in_until(job))

    def omitfrom_jobs(self):
        'Returns a generator of jobs specified by omitfromjobs'
        return (job for job in self.jobs if self.in_omitfrom(job))

    def downstream_of_omitfrom(self):
        "Returns the downstream of --omit-from rules or files."
        return filter(lambda job: not self.in_omitfrom(job),
                      self.bfs(self.depending, *self.omitfrom_jobs()))

    def delete_omitfrom_jobs(self):
        "Removes jobs downstream of jobs specified by --omit-from."
        if not self.omitrules and not self.omitfiles:
            return
        downstream_jobs = list(self.downstream_of_omitfrom()
                               )  # need to cast as list before deleting jobs
        for job in downstream_jobs:
            self.delete_job(job, recursive=False, add_dependencies=True)

    def set_until_jobs(self):
        "Removes jobs downstream of jobs specified by --omit-from."
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

    def update_ready(self):
        """ Update information whether a job is ready to execute. """
        for job in filter(self.needrun, self.jobs):
            if not self.finished(job) and self._ready(job):
                self._ready_jobs.add(job)

    def update_downstream_size(self):
        for job in self.needrun_jobs:
            self._downstream_size[job] = sum(
                1
                for _ in self.bfs(self.depending,
                                  job,
                                  stop=self.noneedrun_finished)) - 1

    def update_temp_input_count(self):
        for job in self.needrun_jobs:
            self._temp_input_count[job] = sum(1 for _ in self.temp_input(job))

    def postprocess(self):
        self.update_needrun()
        self.update_priority()
        self.update_ready()
        self.update_downstream_size()
        self.update_temp_input_count()

    def _ready(self, job):
        return self._finished.issuperset(filter(self.needrun,
                                                self.dependencies[job]))

    def finish(self, job, update_dynamic=True):
        self._finished.add(job)
        try:
            self._ready_jobs.remove(job)
        except KeyError:
            pass
        # mark depending jobs as ready
        for job_ in self.depending[job]:
            if self.needrun(job_) and self._ready(job_):
                self._ready_jobs.add(job_)

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

    def update_dynamic(self, job):
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
        newjob = Job(newrule, self, format_wildcards=non_dynamic_wildcards)
        self.replace_job(job, newjob)
        for job_ in depending:
            if job_.dynamic_input:
                newrule_ = job_.rule.dynamic_branch(dynamic_wildcards)
                if newrule_ is not None:
                    self.specialize_rule(job_.rule, newrule_)
                    if not self.dynamic(job_):
                        logger.debug("Updating job {}.".format(job_))
                        newjob_ = Job(newrule_,
                                      self,
                                      targetfile=job_.targetfile)

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
        if job in self.targetjobs:
            self.targetjobs.remove(job)
            self.targetjobs.add(newjob)
        depending = list(self.depending[job].items())
        if self.finished(job):
            self._finished.add(newjob)

        self.delete_job(job)
        self.update([newjob])

        for job_, files in depending:
            if not job_.dynamic_input:
                self.dependencies[job_][newjob].update(files)
                self.depending[newjob][job_].update(files)

    def specialize_rule(self, rule, newrule):
        assert newrule is not None
        self.rules.add(newrule)
        self.update_output_index()

    def collect_potential_dependencies(self, job):
        dependencies = defaultdict(list)
        # use a set to circumvent multiple jobs for the same file
        # if user specified it twice
        file2jobs = self.file2jobs
        for file in set(job.input):
            # omit the file if it comes from a subworkflow
            if file in job.subworkflow_input:
                continue
            try:
                if file in job.dependencies:
                    jobs = [Job(job.dependencies[file], self, targetfile=file)]
                else:
                    jobs = file2jobs(file)
                dependencies[file].extend(jobs)
            except MissingRuleException as ex:
                pass
        return dependencies

    def bfs(self, direction, *jobs, stop=lambda job: False):
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
        visited = set()
        for job in jobs:
            for job_ in self._dfs(direction,
                                  job,
                                  visited,
                                  stop=stop,
                                  post=post):
                yield job_

    def _dfs(self, direction, job, visited, stop, post):
        if stop(job):
            return
        if not post:
            yield job
        for job_ in direction[job]:
            if not job_ in visited:
                visited.add(job_)
                for j in self._dfs(direction, job_, visited, stop, post):
                    yield j
        if post:
            yield job

    def is_isomorph(self, job1, job2):
        if job1.rule != job2.rule:
            return False
        rule = lambda job: job.rule.name
        queue1, queue2 = [job1], [job2]
        visited1, visited2 = set(queue1), set(queue2)
        while queue1 and queue2:
            job1, job2 = queue1.pop(0), queue2.pop(0)
            deps1 = sorted(self.dependencies[job1], key=rule)
            deps2 = sorted(self.dependencies[job2], key=rule)
            for job1_, job2_ in zip(deps1, deps2):
                if job1_.rule != job2_.rule:
                    return False
                if not job1_ in visited1 and not job2_ in visited2:
                    queue1.append(job1_)
                    visited1.add(job1_)
                    queue2.append(job2_)
                    visited2.add(job2_)
                elif not (job1_ in visited1 and job2_ in visited2):
                    return False
        return True

    def all_longest_paths(self, *jobs):
        paths = defaultdict(list)

        def all_longest_paths(_jobs):
            for job in _jobs:
                if job in paths:
                    continue
                deps = self.dependencies[job]
                if not deps:
                    paths[job].append([job])
                    continue
                all_longest_paths(deps)
                for _job in deps:
                    paths[job].extend(path + [job] for path in paths[_job])

        all_longest_paths(jobs)
        return chain(*(paths[job] for job in jobs))

    def new_wildcards(self, job):
        new_wildcards = set(job.wildcards.items())
        for job_ in self.dependencies[job]:
            if not new_wildcards:
                return set()
            for wildcard in job_.wildcards.items():
                new_wildcards.discard(wildcard)
        return new_wildcards

    def rule2job(self, targetrule):
        return Job(targetrule, self)

    def file2jobs(self, targetfile):
        rules = self.output_index.match(targetfile)
        jobs = []
        exceptions = list()
        for rule in rules:
            if rule.is_producer(targetfile):
                try:
                    jobs.append(Job(rule, self, targetfile=targetfile))
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
            yield "output_file\tdate\trule\tversion\tinput_file(s)\tshellcmd\tstatus\tplan"
        else:
            yield "output_file\tdate\trule\tversion\tstatus\tplan"

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
                    yield "\t".join((f, date, rule, version, input, shellcmd,
                                     status, pending))
                else:
                    yield "\t".join((f, date, rule, version, status, pending))

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
