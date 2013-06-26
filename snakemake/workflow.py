# -*- coding: utf-8 -*-

__author__ = "Johannes Köster"

import re
import os
import sys
import signal
from collections import OrderedDict
from itertools import filterfalse, chain
from operator import attrgetter

from snakemake.logging import logger
from snakemake.rules import Rule, Ruleorder
from snakemake.exceptions import RuleException, CreateRuleException, \
    UnknownRuleException, NoRulesException, print_exception
from snakemake.shell import shell
from snakemake.dag import DAG
from snakemake.scheduler import JobScheduler
from snakemake.parser import parse
from snakemake.io import protected, temp, temporary, expand, dynamic
from snakemake.persistence import Persistence


class Workflow:
    def __init__(self, snakefile=None, snakemakepath=None, jobscript=None):
        """
        Create the controller.
        """
        self._rules = OrderedDict()
        self.first_rule = None
        self._workdir = None
        self._ruleorder = Ruleorder()
        self.linemaps = dict()
        self.rule_count = 0
        self.snakefile = snakefile
        self.snakemakepath = os.path.abspath(snakemakepath)
        self.jobscript = jobscript
        self.persistence = None
        self.global_resources = None
        self.globals = globals()

    @property
    def rules(self):
        return self._rules.values()

    def check(self):
        for clause in self._ruleorder:
            for rulename in clause:
                if not self.is_rule(rulename):
                    raise UnknownRuleException(
                        rulename, prefix = "Error in ruleorder definition.")

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

    def list_rules(self, details=True, log=logger.info):
        log("Available rules:")
        for rule in self.rules:
            log(rule.name)
            if details:
                if rule.docstring:
                    for line in rule.docstring.split("\n"):
                        log("\t" + line)

    def execute(
        self, targets=None, dryrun=False,  touch=False, cores=1,
        forcetargets=False, forceall=False, forcerun=None,
        prioritytargets=None, quiet=False, keepgoing=False,
        printshellcmds=False, printreason=False, printdag=False,
        cluster=None, immediate_submit=False, ignore_ambiguity=False,
        workdir=None, printruledag=False,
        stats=None, force_incomplete=False, ignore_incomplete=False,
        list_version_changes=False, list_code_changes=False,
        list_input_changes=False, list_params_changes=False,
        summary=False, output_wait=3, nolock=False, unlock=False,
        resources=None, notemp=False,
        cleanup_metadata=None):

        self.global_resources = dict() if cluster or resources is None else resources
        self.global_resources["_cores"] = cores

        def rules(items):
            return map(self._rules.__getitem__, filter(self.is_rule, items))

        def files(items):
            return map(os.path.relpath, filterfalse(self.is_rule, items))

        if workdir is None:
            workdir = os.getcwd() if self._workdir is None else self._workdir
        os.chdir(workdir)

        if not targets:
            targets = [self.first_rule]
        if prioritytargets is None:
            prioritytargets = list()
        if forcerun is None:
            forcerun = list()

        priorityrules = set(rules(prioritytargets))
        priorityfiles = set(files(prioritytargets))
        forcerules = set(rules(forcerun))
        forcefiles = set(files(forcerun))
        targetrules = set(chain(
            rules(targets), filterfalse(Rule.has_wildcards, priorityrules),
            filterfalse(Rule.has_wildcards, forcerules)))
        targetfiles = set(chain(files(targets), priorityfiles, forcefiles))
        if forcetargets:
            forcefiles.update(targetfiles)
            forcerules.update(targetrules)

        dag = DAG(
            self, dryrun=dryrun, targetfiles=targetfiles,
            targetrules=targetrules,
            forceall=forceall, forcefiles=forcefiles,
            forcerules=forcerules, priorityfiles=priorityfiles,
            priorityrules=priorityrules, ignore_ambiguity=ignore_ambiguity,
            force_incomplete=force_incomplete,
            ignore_incomplete=ignore_incomplete, notemp=notemp)

        self.persistence = Persistence(nolock=nolock, dag=dag)

        if cleanup_metadata:
            for f in cleanup_metadata:
                self.persistence.cleanup_metadata(f)
            return True

        dag.init()
        dag.check_dynamic()

        if unlock:
            try:
                self.persistence.cleanup_locks()
                logger.warning("Unlocking working directory.")
                return True
            except IOError:
                logger.error("Error: Unlocking the directory {} failed. Maybe "
                "you don't have the permissions?")
                return False
        try:
            self.persistence.lock()
        except IOError:
            logger.critical("Error: Directory cannot be locked. Please make "
                "sure that no other Snakemake process is trying to create "
                "the same files in the following directory:\n{}\n"
                "If you are sure that no other "
                "instances of snakemake are running on this directory, "
                "the remaining lock was likely caused by a kill signal or "
                "a power loss. It can be removed with "
                "the --unlock argument.".format(os.getcwd()))
            return False

        dag.check_incomplete()
        dag.postprocess()

        if printdag:
            print(dag)
            return True
        elif printruledag:
            print(dag.rule_dot())
            return True
        elif summary:
            print("\n".join(dag.summary()))
            return True
        elif list_version_changes:
            items = list(chain(
                *map(self.persistence.version_changed, dag.jobs)))
            if items:
                print(*items, sep="\n")
            return True
        elif list_code_changes:
            items = list(chain(
                *map(self.persistence.code_changed, dag.jobs)))
            if items:
                print(*items, sep="\n")
            return True
        elif list_input_changes:
            items = list(chain(
                *map(self.persistence.input_changed, dag.jobs)))
            if items:
                print(*items, sep="\n")
            return True
        elif list_params_changes:
            items = list(chain(
                *map(self.persistence.params_changed, dag.jobs)))
            if items:
                print(*items, sep="\n")
            return True

        scheduler = JobScheduler(
            self, dag, cores, dryrun=dryrun, touch=touch, cluster=cluster,
            immediate_submit=immediate_submit,
            quiet=quiet, keepgoing=keepgoing,
            printreason=printreason, printshellcmds=printshellcmds,
            output_wait=output_wait)

        if not dryrun and not quiet and len(dag):
            if cluster:
                logger.warning("Provided cluster nodes: {}".format(cores))
            else:
                logger.warning("Provided cores: {}".format(cores))
            logger.warning("\n".join(dag.stats()))

        success = scheduler.schedule()

        if success:
            if dryrun:
                if not quiet:
                    logger.warning("\n".join(dag.stats()))
            elif stats:
                scheduler.stats.to_csv(stats)
        else:
            logger.critical(
                "Exiting because a job execution failed. "
                "Look above for error message")
            return False
        return True

    def include(self, snakefile, workdir=None, overwrite_first_rule=False,
        print_compilation=False):
        """
        Include a snakefile.
        """
        global workflow
        workflow = self
        first_rule = self.first_rule
        if workdir:
            os.chdir(workdir)
        code, linemap = parse(snakefile)

        if print_compilation:
            print(code)

        self.linemaps[snakefile] = linemap
        exec(compile(code, snakefile, "exec"), self.globals)
        if not overwrite_first_rule:
            self.first_rule = first_rule

    def workdir(self, workdir):
        if self._workdir is None:
            if not os.path.exists(workdir):
                os.makedirs(workdir)
            self._workdir = workdir

    def ruleorder(self, *rulenames):
        self._ruleorder.add(*rulenames)

    def rule(self, name=None, lineno=None, snakefile=None):
        name = self.add_rule(name, lineno, snakefile)
        rule = self.get_rule(name)

        def decorate(ruleinfo):
            if ruleinfo.input:
                rule.set_input(*ruleinfo.input[0], **ruleinfo.input[1])
            if ruleinfo.output:
                rule.set_output(*ruleinfo.output[0], **ruleinfo.output[1])
            if ruleinfo.params:
                rule.set_params(*ruleinfo.params[0], **ruleinfo.params[1])
            if ruleinfo.threads:
                if not isinstance(ruleinfo.threads, int):
                    raise RuleException("Threads value has to be an integer.",
                        rule=rule)
                rule.resources["_cores"] = ruleinfo.threads
            if ruleinfo.resources:
                args, resources = ruleinfo.resources
                if args:
                    raise RuleException("Resources have to be named.")
                if not all(map(lambda r: isinstance(r, int), resources.values())):
                    raise RuleException("Resources values have to be integers.", rule=rule)
                rule.resources.update(resources)
            if ruleinfo.priority:
                if (not isinstance(ruleinfo.priority, int)
                    and not isinstance(ruleinfo.priority, float)):
                    raise RuleException("Priority values have to be numeric.",
                        rule=rule)
                rule.priority = ruleinfo.priority
            if ruleinfo.version:
                rule.version = ruleinfo.version
            if ruleinfo.log:
                rule.log = ruleinfo.log
            if ruleinfo.message:
                rule.message = ruleinfo.message
            rule.docstring = ruleinfo.docstring
            rule.run_func = ruleinfo.func
            rule.shellcmd = ruleinfo.shellcmd
            ruleinfo.func.__name__ = "__{}".format(name)
            self.globals[ruleinfo.func.__name__] = ruleinfo.func
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

    def message(self, message):
        def decorate(ruleinfo):
            ruleinfo.message = message
            return ruleinfo
        return decorate

    def threads(self, threads):
        def decorate(ruleinfo):
            ruleinfo.threads = threads
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

    def log(self, log):
        def decorate(ruleinfo):
            ruleinfo.log = log
            return ruleinfo
        return decorate

    def shellcmd(self, cmd):
        def decorate(ruleinfo):
            ruleinfo.shellcmd = cmd
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
        self.input = None
        self.output = None
        self.params = None
        self.message = None
        self.threads = None
        self.resources = None
        self.priority = None
        self.version = None
        self.log = None
        self.docstring = None
