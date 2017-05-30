__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import re
import os
import sys
import signal
import json
import urllib
from collections import OrderedDict
from itertools import filterfalse, chain
from functools import partial
from operator import attrgetter
import copy

from snakemake.logging import logger, format_resources, format_resource_names
from snakemake.rules import Rule, Ruleorder
from snakemake.exceptions import RuleException, CreateRuleException, \
    UnknownRuleException, NoRulesException, print_exception, WorkflowError
from snakemake.shell import shell
from snakemake.dag import DAG
from snakemake.scheduler import JobScheduler
from snakemake.parser import parse
import snakemake.io
from snakemake.io import protected, temp, temporary, ancient, expand, dynamic, glob_wildcards, flag, not_iterable, touch, unpack
from snakemake.persistence import Persistence
from snakemake.utils import update_config
from snakemake.script import script
from snakemake.wrapper import wrapper
import snakemake.wrapper
from snakemake.common import Mode

class Workflow:
    def __init__(self,
                 snakefile=None,
                 jobscript=None,
                 overwrite_shellcmd=None,
                 overwrite_config=dict(),
                 overwrite_workdir=None,
                 overwrite_configfile=None,
                 overwrite_clusterconfig=dict(),
                 config_args=None,
                 debug=False,
                 use_conda=False,
                 conda_prefix=None,
                 mode=Mode.default,
                 wrapper_prefix=None,
                 printshellcmds=False,
                 restart_times=None):
        """
        Create the controller.
        """
        self._rules = OrderedDict()
        self.first_rule = None
        self._workdir = None
        self.overwrite_workdir = overwrite_workdir
        self.workdir_init = os.path.abspath(os.curdir)
        self._ruleorder = Ruleorder()
        self._localrules = set()
        self.linemaps = dict()
        self.rule_count = 0
        self.basedir = os.path.dirname(snakefile)
        self.snakefile = os.path.abspath(snakefile)
        self.included = []
        self.included_stack = []
        self.jobscript = jobscript
        self.persistence = None
        self.global_resources = None
        self.globals = globals()
        self._subworkflows = dict()
        self.overwrite_shellcmd = overwrite_shellcmd
        self.overwrite_config = overwrite_config
        self.overwrite_configfile = overwrite_configfile
        self.overwrite_clusterconfig = overwrite_clusterconfig
        self.config_args = config_args
        self.immediate_submit = None
        self._onsuccess = lambda log: None
        self._onerror = lambda log: None
        self._onstart = lambda log: None
        self._wildcard_constraints = dict()
        self.debug = debug
        self._rulecount = 0
        self.use_conda = use_conda
        self.conda_prefix = conda_prefix
        self.mode = mode
        self.wrapper_prefix = wrapper_prefix
        self.printshellcmds = printshellcmds
        self.restart_times = restart_times

        global config
        config = copy.deepcopy(self.overwrite_config)

        global cluster_config
        cluster_config = copy.deepcopy(self.overwrite_clusterconfig)

        global rules
        rules = Rules()

    @property
    def subworkflows(self):
        return self._subworkflows.values()

    @property
    def rules(self):
        return self._rules.values()

    @property
    def concrete_files(self):
        return (
            file
            for rule in self.rules for file in chain(rule.input, rule.output)
            if not callable(file) and not file.contains_wildcard()
        )

    def check(self):
        for clause in self._ruleorder:
            for rulename in clause:
                if not self.is_rule(rulename):
                    raise UnknownRuleException(
                        rulename,
                        prefix="Error in ruleorder definition.")

    def add_rule(self, name=None, lineno=None, snakefile=None):
        """
        Add a rule.
        """
        if name is None:
            name = str(len(self._rules) + 1)
        if self.is_rule(name):
            raise CreateRuleException(
                "The name {} is already used by another rule".format(name))
        rule = Rule(name, self, lineno=lineno, snakefile=snakefile)
        self._rules[rule.name] = rule
        self.rule_count += 1
        if not self.first_rule:
            self.first_rule = rule.name
        return name

    def is_rule(self, name):
        """
        Return True if name is the name of a rule.

        Arguments
        name -- a name
        """
        return name in self._rules

    def get_rule(self, name):
        """
        Get rule by name.

        Arguments
        name -- the name of the rule
        """
        if not self._rules:
            raise NoRulesException()
        if not name in self._rules:
            raise UnknownRuleException(name)
        return self._rules[name]

    def list_rules(self, only_targets=False):
        rules = self.rules
        if only_targets:
            rules = filterfalse(Rule.has_wildcards, rules)
        for rule in rules:
            logger.rule_info(name=rule.name, docstring=rule.docstring)

    def list_resources(self):
        for resource in set(
            resource for rule in self.rules for resource in rule.resources):
            if resource not in "_cores _nodes".split():
                logger.info(resource)

    def is_local(self, rule):
        return rule.name in self._localrules or rule.norun

    def execute(self,
                targets=None,
                dryrun=False,
                touch=False,
                cores=1,
                nodes=1,
                local_cores=1,
                forcetargets=False,
                forceall=False,
                forcerun=None,
                until=[],
                omit_from=[],
                prioritytargets=None,
                quiet=False,
                keepgoing=False,
                printshellcmds=False,
                printreason=False,
                printdag=False,
                cluster=None,
                cluster_sync=None,
                jobname=None,
                immediate_submit=False,
                ignore_ambiguity=False,
                printrulegraph=False,
                printd3dag=False,
                drmaa=None,
                drmaa_log_dir=None,
                stats=None,
                force_incomplete=False,
                ignore_incomplete=False,
                list_version_changes=False,
                list_code_changes=False,
                list_input_changes=False,
                list_params_changes=False,
                summary=False,
                archive=None,
                detailed_summary=False,
                latency_wait=3,
                benchmark_repeats=3,
                wait_for_files=None,
                nolock=False,
                unlock=False,
                resources=None,
                notemp=False,
                nodeps=False,
                cleanup_metadata=None,
                subsnakemake=None,
                updated_files=None,
                keep_target_files=False,
                keep_shadow=False,
                keep_remote_local=False,
                allowed_rules=None,
                max_jobs_per_second=None,
                greediness=1.0,
                no_hooks=False,
                force_use_threads=False):

        self.global_resources = dict() if resources is None else resources
        self.global_resources["_cores"] = cores
        self.global_resources["_nodes"] = nodes
        self.immediate_submit = immediate_submit

        def rules(items):
            return map(self._rules.__getitem__, filter(self.is_rule, items))

        if keep_target_files:

            def files(items):
                return filterfalse(self.is_rule, items)
        else:

            def files(items):
                return map(os.path.relpath, filterfalse(self.is_rule, items))

        if not targets:
            targets = [self.first_rule
                       ] if self.first_rule is not None else list()

        if prioritytargets is None:
            prioritytargets = list()
        if forcerun is None:
            forcerun = list()
        if until is None:
            until = list()
        if omit_from is None:
            omit_from = list()

        priorityrules = set(rules(prioritytargets))
        priorityfiles = set(files(prioritytargets))
        forcerules = set(rules(forcerun))
        forcefiles = set(files(forcerun))
        untilrules = set(rules(until))
        untilfiles = set(files(until))
        omitrules = set(rules(omit_from))
        omitfiles = set(files(omit_from))
        targetrules = set(chain(rules(targets),
                                filterfalse(Rule.has_wildcards, priorityrules),
                                filterfalse(Rule.has_wildcards, forcerules),
                                filterfalse(Rule.has_wildcards, untilrules)))
        targetfiles = set(chain(files(targets), priorityfiles, forcefiles, untilfiles))
        if forcetargets:
            forcefiles.update(targetfiles)
            forcerules.update(targetrules)

        rules = self.rules
        if allowed_rules:
            rules = [rule for rule in rules if rule.name in set(allowed_rules)]

        if wait_for_files is not None:
            try:
                snakemake.io.wait_for_files(wait_for_files,
                                            latency_wait=latency_wait)
            except IOError as e:
                logger.error(str(e))
                return False

        dag = DAG(
            self, rules,
            dryrun=dryrun,
            targetfiles=targetfiles,
            targetrules=targetrules,
            forceall=forceall,
            forcefiles=forcefiles,
            forcerules=forcerules,
            priorityfiles=priorityfiles,
            priorityrules=priorityrules,
            untilfiles=untilfiles,
            untilrules=untilrules,
            omitfiles=omitfiles,
            omitrules=omitrules,
            ignore_ambiguity=ignore_ambiguity,
            force_incomplete=force_incomplete,
            ignore_incomplete=ignore_incomplete or printdag or printrulegraph,
            notemp=notemp,
            keep_remote_local=keep_remote_local)

        self.persistence = Persistence(
            nolock=nolock,
            dag=dag,
            conda_prefix=self.conda_prefix,
            warn_only=dryrun or printrulegraph or printdag or summary or archive or
            list_version_changes or list_code_changes or list_input_changes or
            list_params_changes)

        if cleanup_metadata:
            for f in cleanup_metadata:
                self.persistence.cleanup_metadata(f)
            return True

        dag.init()
        dag.check_dynamic()

        if unlock:
            try:
                self.persistence.cleanup_locks()
                logger.info("Unlocking working directory.")
                return True
            except IOError:
                logger.error("Error: Unlocking the directory {} failed. Maybe "
                             "you don't have the permissions?")
                return False
        try:
            self.persistence.lock()
        except IOError:
            logger.error(
                "Error: Directory cannot be locked. Please make "
                "sure that no other Snakemake process is trying to create "
                "the same files in the following directory:\n{}\n"
                "If you are sure that no other "
                "instances of snakemake are running on this directory, "
                "the remaining lock was likely caused by a kill signal or "
                "a power loss. It can be removed with "
                "the --unlock argument.".format(os.getcwd()))
            return False

        if self.subworkflows and not printdag and not printrulegraph:
            # backup globals
            globals_backup = dict(self.globals)
            # execute subworkflows
            for subworkflow in self.subworkflows:
                subworkflow_targets = subworkflow.targets(dag)
                updated = list()
                if subworkflow_targets:
                    logger.info(
                        "Executing subworkflow {}.".format(subworkflow.name))
                    if not subsnakemake(subworkflow.snakefile,
                                        workdir=subworkflow.workdir,
                                        targets=subworkflow_targets,
                                        configfile=subworkflow.configfile,
                                        updated_files=updated):
                        return False
                    dag.updated_subworkflow_files.update(subworkflow.target(f)
                                                         for f in updated)
                else:
                    logger.info("Subworkflow {}: Nothing to be done.".format(
                        subworkflow.name))
            if self.subworkflows:
                logger.info("Executing main workflow.")
            # rescue globals
            self.globals.update(globals_backup)

        dag.check_incomplete()
        dag.postprocess()

        if nodeps:
            missing_input = [f for job in dag.targetjobs for f in job.input
                             if dag.needrun(job) and not os.path.exists(f)]
            if missing_input:
                logger.error(
                    "Dependency resolution disabled (--nodeps) "
                    "but missing input "
                    "files detected. If this happens on a cluster, please make sure "
                    "that you handle the dependencies yourself or turn of "
                    "--immediate-submit. Missing input files:\n{}".format(
                        "\n".join(missing_input)))
                return False

        updated_files.extend(f for job in dag.needrun_jobs for f in job.output)

        if printd3dag:
            dag.d3dag()
            return True
        elif printdag:
            print(dag)
            return True
        elif printrulegraph:
            print(dag.rule_dot())
            return True
        elif summary:
            print("\n".join(dag.summary(detailed=False)))
            return True
        elif detailed_summary:
            print("\n".join(dag.summary(detailed=True)))
            return True
        elif archive:
            dag.archive(archive)
            return True
        elif list_version_changes:
            items = list(
                chain(*map(self.persistence.version_changed, dag.jobs)))
            if items:
                print(*items, sep="\n")
            return True
        elif list_code_changes:
            items = list(chain(*map(self.persistence.code_changed, dag.jobs)))
            for j in dag.jobs:
                items.extend(list(j.outputs_older_than_script()))
            if items:
                print(*items, sep="\n")
            return True
        elif list_input_changes:
            items = list(chain(*map(self.persistence.input_changed, dag.jobs)))
            if items:
                print(*items, sep="\n")
            return True
        elif list_params_changes:
            items = list(
                chain(*map(self.persistence.params_changed, dag.jobs)))
            if items:
                print(*items, sep="\n")
            return True

        if not keep_shadow:
            self.persistence.cleanup_shadow()

        if self.use_conda:
            dag.create_conda_envs(dryrun=dryrun)

        scheduler = JobScheduler(self, dag, cores,
                                 local_cores=local_cores,
                                 dryrun=dryrun,
                                 touch=touch,
                                 cluster=cluster,
                                 cluster_config=cluster_config,
                                 cluster_sync=cluster_sync,
                                 jobname=jobname,
                                 max_jobs_per_second=max_jobs_per_second,
                                 quiet=quiet,
                                 keepgoing=keepgoing,
                                 drmaa=drmaa,
                                 drmaa_log_dir=drmaa_log_dir,
                                 printreason=printreason,
                                 printshellcmds=printshellcmds,
                                 latency_wait=latency_wait,
                                 benchmark_repeats=benchmark_repeats,
                                 greediness=greediness,
                                 force_use_threads=force_use_threads)

        if not dryrun:
            if len(dag):
                if cluster or cluster_sync or drmaa:
                    logger.resources_info(
                        "Provided cluster nodes: {}".format(nodes))
                else:
                    logger.resources_info("Provided cores: {}".format(cores))
                    logger.resources_info("Rules claiming more threads will be scaled down.")
                provided_resources = format_resources(resources)
                if provided_resources:
                    logger.resources_info(
                        "Provided resources: " + provided_resources)
                unlimited_resources = format_resource_names(set(
                    resource for job in dag.needrun_jobs
                    for resource in job.resources.keys()
                    if resource not in resources))
                if unlimited_resources:
                    logger.resources_info(
                        "Unlimited resources: " + unlimited_resources)
                logger.run_info("\n".join(dag.stats()))
            else:
                logger.info("Nothing to be done.")
        if dryrun and not len(dag):
            logger.info("Nothing to be done.")

        if not dryrun and not no_hooks:
            self._onstart(logger.get_logfile())

        success = scheduler.schedule()

        if success:
            if dryrun:
                if len(dag):
                    logger.run_info("\n".join(dag.stats()))
            elif stats:
                scheduler.stats.to_json(stats)
            if not dryrun and not no_hooks:
                self._onsuccess(logger.get_logfile())
            return True
        else:
            if not dryrun and not no_hooks:
                self._onerror(logger.get_logfile())
            return False

    @property
    def current_basedir(self):
        """Basedir of currently parsed Snakefile."""
        assert self.included_stack
        return os.path.abspath(os.path.dirname(self.included_stack[-1]))

    def include(self, snakefile,
                overwrite_first_rule=False,
                print_compilation=False,
                overwrite_shellcmd=None):
        """
        Include a snakefile.
        """
        # check if snakefile is a path to the filesystem
        if not urllib.parse.urlparse(snakefile).scheme:
            if not os.path.isabs(snakefile) and self.included_stack:
                snakefile = os.path.join(self.current_basedir, snakefile)
            # Could still be an url if relative import was used
            if not urllib.parse.urlparse(snakefile).scheme:
                snakefile = os.path.abspath(snakefile)
        # else it could be an url.
        # at least we don't want to modify the path for clarity.

        if snakefile in self.included:
            logger.info("Multiple include of {} ignored".format(snakefile))
            return
        self.included.append(snakefile)
        self.included_stack.append(snakefile)

        global workflow

        workflow = self

        first_rule = self.first_rule
        code, linemap, rulecount = parse(snakefile,
                                         overwrite_shellcmd=self.overwrite_shellcmd,
                                         rulecount=self._rulecount)
        self._rulecount = rulecount

        if print_compilation:
            print(code)

        # insert the current directory into sys.path
        # this allows to import modules from the workflow directory
        sys.path.insert(0, os.path.dirname(snakefile))

        self.linemaps[snakefile] = linemap
        exec(compile(code, snakefile, "exec"), self.globals)
        if not overwrite_first_rule:
            self.first_rule = first_rule
        self.included_stack.pop()

    def onstart(self, func):
        """Register onstart function."""
        self._onstart = func

    def onsuccess(self, func):
        """Register onsuccess function."""
        self._onsuccess = func

    def onerror(self, func):
        """Register onerror function."""
        self._onerror = func

    def global_wildcard_constraints(self, **content):
        """Register global wildcard constraints."""
        self._wildcard_constraints.update(content)

    def workdir(self, workdir):
        """Register workdir."""
        if self.overwrite_workdir is None:
            os.makedirs(workdir, exist_ok=True)
            self._workdir = workdir
            os.chdir(workdir)

    def configfile(self, jsonpath):
        """ Update the global config with the given dictionary. """
        global config
        c = snakemake.io.load_configfile(jsonpath)
        update_config(config, c)
        update_config(config, self.overwrite_config)

    def ruleorder(self, *rulenames):
        self._ruleorder.add(*rulenames)

    def subworkflow(self, name, snakefile=None, workdir=None, configfile=None):
        # Take absolute path of config file, because it is relative to current
        # workdir, which could be changed for the subworkflow.
        if configfile:
            configfile = os.path.abspath(configfile)
        sw = Subworkflow(self, name, snakefile, workdir, configfile)
        self._subworkflows[name] = sw
        self.globals[name] = sw.target

    def localrules(self, *rulenames):
        self._localrules.update(rulenames)

    def rule(self, name=None, lineno=None, snakefile=None):
        name = self.add_rule(name, lineno, snakefile)
        rule = self.get_rule(name)

        def decorate(ruleinfo):
            if ruleinfo.wildcard_constraints:
                rule.set_wildcard_constraints(*ruleinfo.wildcard_constraints[0], **ruleinfo.wildcard_constraints[1])
            if ruleinfo.input:
                rule.set_input(*ruleinfo.input[0], **ruleinfo.input[1])
            if ruleinfo.output:
                rule.set_output(*ruleinfo.output[0], **ruleinfo.output[1])
            if ruleinfo.params:
                rule.set_params(*ruleinfo.params[0], **ruleinfo.params[1])
            if ruleinfo.threads:
                if not isinstance(ruleinfo.threads, int) and not callable(ruleinfo.threads):
                    raise RuleException("Threads value has to be an integer or a callable.",
                                        rule=rule)
                rule.resources["_cores"] = ruleinfo.threads
            if ruleinfo.shadow_depth:
                if ruleinfo.shadow_depth not in (True, "shallow", "full"):
                    raise RuleException(
                        "Shadow must either be 'shallow', 'full', "
                        "or True (equivalent to 'full')", rule=rule)
                if ruleinfo.shadow_depth is True:
                    rule.shadow_depth = 'full'
                else:
                    rule.shadow_depth = ruleinfo.shadow_depth
            if ruleinfo.resources:
                args, resources = ruleinfo.resources
                if args:
                    raise RuleException("Resources have to be named.")
                if not all(map(lambda r: isinstance(r, int) or callable(r),
                               resources.values())):
                    raise RuleException(
                        "Resources values have to be integers or callables",
                        rule=rule)
                rule.resources.update(resources)
            if ruleinfo.priority:
                if (not isinstance(ruleinfo.priority, int) and
                    not isinstance(ruleinfo.priority, float)):
                    raise RuleException("Priority values have to be numeric.",
                                        rule=rule)
                rule.priority = ruleinfo.priority
            if ruleinfo.version:
                rule.version = ruleinfo.version
            if ruleinfo.log:
                rule.set_log(*ruleinfo.log[0], **ruleinfo.log[1])
            if ruleinfo.message:
                rule.message = ruleinfo.message
            if ruleinfo.benchmark:
                rule.benchmark = ruleinfo.benchmark
            if ruleinfo.wrapper:
                rule.conda_env = snakemake.wrapper.get_conda_env(ruleinfo.wrapper)
            if ruleinfo.conda_env:
                if not (ruleinfo.script or ruleinfo.wrapper or ruleinfo.shellcmd):
                    raise RuleException("Conda environments are only allowed "
                        "with shell, script or wrapper directives (not with run).", rule=rule)
                if not os.path.isabs(ruleinfo.conda_env):
                    ruleinfo.conda_env = os.path.join(self.current_basedir, ruleinfo.conda_env)
                rule.conda_env = ruleinfo.conda_env
            rule.norun = ruleinfo.norun
            rule.docstring = ruleinfo.docstring
            rule.run_func = ruleinfo.func
            rule.shellcmd = ruleinfo.shellcmd
            rule.script = ruleinfo.script
            rule.wrapper = ruleinfo.wrapper
            rule.restart_times=self.restart_times

            ruleinfo.func.__name__ = "__{}".format(name)
            self.globals[ruleinfo.func.__name__] = ruleinfo.func
            setattr(rules, name, rule)
            return ruleinfo.func

        return decorate

    def docstring(self, string):
        def decorate(ruleinfo):
            ruleinfo.docstring = string
            return ruleinfo

        return decorate

    def input(self, *paths, **kwpaths):
        def decorate(ruleinfo):
            ruleinfo.input = (paths, kwpaths)
            return ruleinfo

        return decorate

    def output(self, *paths, **kwpaths):
        def decorate(ruleinfo):
            ruleinfo.output = (paths, kwpaths)
            return ruleinfo

        return decorate

    def params(self, *params, **kwparams):
        def decorate(ruleinfo):
            ruleinfo.params = (params, kwparams)
            return ruleinfo

        return decorate

    def wildcard_constraints(self, *wildcard_constraints, **kwwildcard_constraints):
        def decorate(ruleinfo):
            ruleinfo.wildcard_constraints = (wildcard_constraints, kwwildcard_constraints)
            return ruleinfo

        return decorate

    def message(self, message):
        def decorate(ruleinfo):
            ruleinfo.message = message
            return ruleinfo

        return decorate

    def benchmark(self, benchmark):
        def decorate(ruleinfo):
            ruleinfo.benchmark = benchmark
            return ruleinfo

        return decorate

    def conda(self, conda_env):
        def decorate(ruleinfo):
            ruleinfo.conda_env = conda_env
            return ruleinfo

        return decorate

    def threads(self, threads):
        def decorate(ruleinfo):
            ruleinfo.threads = threads
            return ruleinfo

        return decorate

    def shadow(self, shadow_depth):
        def decorate(ruleinfo):
            ruleinfo.shadow_depth = shadow_depth
            return ruleinfo

        return decorate

    def resources(self, *args, **resources):
        def decorate(ruleinfo):
            ruleinfo.resources = (args, resources)
            return ruleinfo

        return decorate

    def priority(self, priority):
        def decorate(ruleinfo):
            ruleinfo.priority = priority
            return ruleinfo

        return decorate

    def version(self, version):
        def decorate(ruleinfo):
            ruleinfo.version = version
            return ruleinfo

        return decorate

    def log(self, *logs, **kwlogs):
        def decorate(ruleinfo):
            ruleinfo.log = (logs, kwlogs)
            return ruleinfo

        return decorate

    def shellcmd(self, cmd):
        def decorate(ruleinfo):
            ruleinfo.shellcmd = cmd
            return ruleinfo

        return decorate

    def script(self, script):
        def decorate(ruleinfo):
            ruleinfo.script = script
            return ruleinfo

        return decorate

    def wrapper(self, wrapper):
        def decorate(ruleinfo):
            ruleinfo.wrapper = wrapper
            return ruleinfo

        return decorate

    def norun(self):
        def decorate(ruleinfo):
            ruleinfo.norun = True
            return ruleinfo

        return decorate

    def run(self, func):
        return RuleInfo(func)

    @staticmethod
    def _empty_decorator(f):
        return f


class RuleInfo:
    def __init__(self, func):
        self.func = func
        self.shellcmd = None
        self.norun = False
        self.input = None
        self.output = None
        self.params = None
        self.message = None
        self.benchmark = None
        self.conda_env = None
        self.wildcard_constraints = None
        self.threads = None
        self.shadow_depth = None
        self.resources = None
        self.priority = None
        self.version = None
        self.log = None
        self.docstring = None
        self.script = None
        self.wrapper = None


class Subworkflow:
    def __init__(self, workflow, name, snakefile, workdir, configfile):
        self.workflow = workflow
        self.name = name
        self._snakefile = snakefile
        self._workdir = workdir
        self.configfile = configfile

    @property
    def snakefile(self):
        if self._snakefile is None:
            return os.path.abspath(os.path.join(self.workdir, "Snakefile"))
        if not os.path.isabs(self._snakefile):
            return os.path.abspath(os.path.join(self.workflow.basedir,
                                                self._snakefile))
        return self._snakefile

    @property
    def workdir(self):
        workdir = "." if self._workdir is None else self._workdir
        if not os.path.isabs(workdir):
            return os.path.abspath(os.path.join(self.workflow.basedir,
                                                workdir))
        return workdir

    def target(self, paths):
        if not_iterable(paths):
            return flag(os.path.join(self.workdir, paths), "subworkflow", self)
        return [self.target(path) for path in paths]

    def targets(self, dag):
        return [f for job in dag.jobs for f in job.subworkflow_input
                if job.subworkflow_input[f] is self]


class Rules:
    """ A namespace for rules so that they can be accessed via dot notation. """
    pass


def srcdir(path):
    """Return the absolute path, relative to the source directory of the current Snakefile."""
    if not workflow.included_stack:
        return None
    return os.path.join(os.path.dirname(workflow.included_stack[-1]), path)
