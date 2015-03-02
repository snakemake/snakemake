# -*- coding: utf-8 -*-

import os
import sys
import base64
import json

from collections import defaultdict
from itertools import chain
from functools import partial
from operator import attrgetter

from snakemake.io import IOFile, Wildcards, Resources, _IOFile
from snakemake.utils import format, listfiles
from snakemake.exceptions import RuleException, ProtectedOutputException
from snakemake.exceptions import UnexpectedOutputException
from snakemake.logging import logger

__author__ = "Johannes Köster"


def jobfiles(jobs, type):
    return chain(*map(attrgetter(type), jobs))


class Job:
    HIGHEST_PRIORITY = sys.maxsize

    def __init__(self, rule, dag, targetfile=None, format_wildcards=None):
        self.rule = rule
        self.dag = dag
        self.targetfile = targetfile

        self.wildcards_dict = self.rule.get_wildcards(targetfile)
        self.wildcards = Wildcards(fromdict=self.wildcards_dict)
        self._format_wildcards = (self.wildcards
            if format_wildcards is None
            else Wildcards(fromdict=format_wildcards))

        (
            self.input, self.output, self.params,
            self.log, self.benchmark,
            self.ruleio, self.dependencies
        ) = rule.expand_wildcards(self.wildcards_dict)

        self.resources_dict = {
            name: min(self.rule.workflow.global_resources.get(name, res), res)
            for name, res in rule.resources.items()}
        self.threads = self.resources_dict["_cores"]
        self.resources = Resources(fromdict=self.resources_dict)
        self._inputsize = None

        self.dynamic_output, self.dynamic_input = set(), set()
        self.temp_output, self.protected_output = set(), set()
        self.touch_output = set()
        self.subworkflow_input = dict()
        for f in self.output:
            f_ = self.ruleio[f]
            if f_ in self.rule.dynamic_output:
                self.dynamic_output.add(f)
            if f_ in self.rule.temp_output:
                self.temp_output.add(f)
            if f_ in self.rule.protected_output:
                self.protected_output.add(f)
            if f_ in self.rule.touch_output:
                self.touch_output.add(f)
        for f in self.input:
            f_ = self.ruleio[f]
            if f_ in self.rule.dynamic_input:
                self.dynamic_input.add(f)
            if f_ in self.rule.subworkflow_input:
                self.subworkflow_input[f] = self.rule.subworkflow_input[f_]
        self._hash = self.rule.__hash__()
        if True or not self.dynamic_output:
            for o in self.output:
                self._hash ^= o.__hash__()

    @property
    def priority(self):
        return self.dag.priority(self)

    @property
    def b64id(self):
        return base64.b64encode((self.rule.name +
            "".join(self.output)).encode("utf-8")).decode("utf-8")

    @property
    def inputsize(self):
        """
        Return the size of the input files.
        Input files need to be present.
        """
        if self._inputsize is None:
            self._inputsize = sum(map(os.path.getsize, self.input))
        return self._inputsize

    @property
    def message(self):
        """ Return the message for this job. """
        try:
            return (self.format_wildcards(self.rule.message)
                if self.rule.message else None)
        except AttributeError as ex:
            raise RuleException(str(ex), rule=self.rule)
        except KeyError as ex:
            raise RuleException("Unknown variable in message "
                "of shell command: {}".format(str(ex)), rule=self.rule)

    @property
    def shellcmd(self):
        """ Return the shell command. """
        try:
            return (self.format_wildcards(self.rule.shellcmd)
                if self.rule.shellcmd else None)
        except AttributeError as ex:
            raise RuleException(str(ex), rule=self.rule)
        except KeyError as ex:
            raise RuleException("Unknown variable when printing "
                "shell command: {}".format(str(ex)), rule=self.rule)

    @property
    def expanded_output(self):
        """ Iterate over output files while dynamic output is expanded. """
        for f, f_ in zip(self.output, self.rule.output):
            if f in self.dynamic_output:
                expansion = self.expand_dynamic(
                    f_,
                    restriction=self.wildcards,
                    omit_value=_IOFile.dynamic_fill
                )
                if not expansion:
                    yield f_
                for f, _ in expansion:
                    yield IOFile(f, self.rule)
            else:
                yield f

    @property
    def dynamic_wildcards(self):
        """ Return all wildcard values determined from dynamic output. """
        combinations = set()
        for f, f_ in zip(self.output, self.rule.output):
            if f in self.dynamic_output:
                for f, w in self.expand_dynamic(
                    f_,
                    restriction=self.wildcards,
                    omit_value=_IOFile.dynamic_fill
                ):
                    combinations.add(tuple(w.items()))
        wildcards = defaultdict(list)
        for combination in combinations:
            for name, value in combination:
                wildcards[name].append(value)
        return wildcards

    @property
    def missing_input(self):
        """ Return missing input files. """
        # omit file if it comes from a subworkflow
        return set(f for f in self.input if not f.exists and not f in self.subworkflow_input)

    @property
    def output_mintime(self):
        """ Return oldest output file. """
        existing = [f.mtime for f in self.expanded_output if f.exists]
        if self.benchmark and self.benchmark.exists:
            existing.append(self.benchmark.mtime)
        if existing:
            return min(existing)
        return None

    @property
    def input_maxtime(self):
        """ Return newest input file. """
        existing = [f.mtime for f in self.input if f.exists]
        if existing:
            return max(existing)
        return None

    def missing_output(self, requested=None):
        """ Return missing output files. """
        files = set()
        if self.benchmark and (requested is None or self.benchmark in requested):
            if not self.benchmark.exists:
                files.add(self.benchmark)

        for f, f_ in zip(self.output, self.rule.output):
            if requested is None or f in requested:
                if f in self.dynamic_output:
                    if not self.expand_dynamic(
                        f_,
                        restriction=self.wildcards,
                        omit_value=_IOFile.dynamic_fill
                    ):
                        files.add("{} (dynamic)".format(f_))
                elif not f.exists:
                    files.add(f)
        return files

    @property
    def existing_output(self):
        return filter(lambda f: f.exists, self.expanded_output)

    def check_protected_output(self):
        protected = list(filter(lambda f: f.protected, self.expanded_output))
        if protected:
            raise ProtectedOutputException(self.rule, protected)

    def prepare(self):
        """
        Prepare execution of job.
        This includes creation of directories and deletion of previously
        created dynamic files.
        """

        self.check_protected_output()

        unexpected_output = self.dag.reason(self).missing_output.intersection(
            self.existing_output)
        if unexpected_output:
            raise UnexpectedOutputException(self.rule, unexpected_output)

        if self.dynamic_output:
            for f, _ in chain(*map(
                partial(
                    self.expand_dynamic,
                    restriction=self.wildcards,
                    omit_value=_IOFile.dynamic_fill),
                self.rule.dynamic_output)):
                os.remove(f)
        for f, f_ in zip(self.output, self.rule.output):
            f.prepare()
        if self.benchmark:
            self.benchmark.prepare()
        if self.log:
            self.log.prepare()

    def cleanup(self):
        """ Cleanup output files. """
        to_remove = [f for f in self.expanded_output if f.exists]
        if to_remove:
            logger.info(
                "Removing output files of failed job {}"
                " since they might be corrupted:\n{}".format(
                    self,
                    ", ".join(to_remove)
                )
            )
            for f in to_remove:
                f.remove()

    def format_wildcards(self, string, **variables):
        """ Format a string with variables from the job. """
        _variables = dict()
        _variables.update(self.rule.workflow.globals)
        _variables.update(dict(
            input=self.input,
            output=self.output,
            params=self.params,
            wildcards=self._format_wildcards,
            threads=self.threads,
            resources=self.resources,
            log=self.log,
            version=self.rule.version
        ))
        _variables.update(variables)
        try:
            return format(string, **_variables)
        except NameError as ex:
            raise RuleException("NameError: " + str(ex), rule=self.rule)
        except IndexError as ex:
            raise RuleException("IndexError: " + str(ex), rule=self.rule)

    def properties(self, omit_resources="_cores _nodes".split()):
        resources = {name: res for name, res in self.resources.items() if name not in omit_resources}
        params = {name: value for name, value in self.params.items()}
        properties = {
            "rule": self.rule.name,
            "local": self.dag.workflow.is_local(self.rule),
            "input": self.input,
            "output": self.output,
            "params": params,
            "threads": self.threads,
            "resources": resources
        }
        return properties

    def json(self):
        return json.dumps(self.properties())

    def __repr__(self):
        return self.rule.name

    def __eq__(self, other):
        if other is None:
            return False
        return self.rule == other.rule and (self.dynamic_output
            or self.wildcards_dict == other.wildcards_dict)

    def __lt__(self, other):
        return self.rule.__lt__(other.rule)

    def __gt__(self, other):
        return self.rule.__gt__(other.rule)

    def __hash__(self):
        return self._hash

    @staticmethod
    def expand_dynamic(pattern, restriction=None, omit_value=None):
        """ Expand dynamic files. """
        return list(listfiles(
            pattern, restriction=restriction, omit_value=omit_value))


class Reason:
    def __init__(self):
        self.updated_input = set()
        self.updated_input_run = set()
        self.missing_output = set()
        self.incomplete_output = set()
        self.forced = False
        self.noio = False
        self.nooutput = False
        self.derived = True

    def __str__(self):
        s = list()
        if self.forced:
            s.append("Forced execution")
        else:
            if self.noio:
                s.append("Rules with neither input nor "
                    "output files are always executed.")
            elif self.nooutput:
                s.append("Rules with a run or shell declaration but no output "
                    "are always executed.")
            else:
                if self.missing_output:
                    s.append("Missing output files: {}".format(
                        ", ".join(self.missing_output)))
                if self.incomplete_output:
                    s.append("Incomplete output files: {}".format(
                        ", ".join(self.incomplete_output)))
                updated_input = self.updated_input - self.updated_input_run
                if updated_input:
                    s.append("Updated input files: {}".format(
                        ", ".join(updated_input)))
                if self.updated_input_run:
                    s.append("This run updates input files: {}".format(
                        ", ".join(self.updated_input_run)))
        s = "; ".join(s)
        #if not self.derived:
        #    s += " (root)"
        return s

    def __bool__(self):
        return bool(self.updated_input or self.missing_output or self.forced
            or self.updated_input_run or self.noio or self.nooutput)
