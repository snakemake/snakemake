__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import hashlib
import os
import sys
import base64
import tempfile
import subprocess

from collections import defaultdict
from itertools import chain
from functools import partial
from operator import attrgetter
from urllib.request import urlopen
from urllib.parse import urlparse

from snakemake.io import IOFile, Wildcards, Resources, _IOFile, is_flagged, contains_wildcard, lstat
from snakemake.utils import format, listfiles
from snakemake.exceptions import RuleException, ProtectedOutputException, WorkflowError
from snakemake.exceptions import UnexpectedOutputException, CreateCondaEnvironmentException
from snakemake.logging import logger
from snakemake.common import DYNAMIC_FILL, lazy_property
from snakemake import conda, wrapper


def jobfiles(jobs, type):
    return chain(*map(attrgetter(type), jobs))


class Job:
    HIGHEST_PRIORITY = sys.maxsize

    __slots__ = ["rule", "dag", "wildcards_dict", "wildcards",
                 "_format_wildcards", "input", "dependencies", "output",
                 "_params", "_log", "_benchmark", "_resources",
                 "_conda_env_file", "_conda_env", "shadow_dir", "_inputsize",
                 "restart_times", "dynamic_output", "dynamic_input",
                 "temp_output", "protected_output", "touch_output",
                 "subworkflow_input", "_hash"]

    def __init__(self, rule, dag, wildcards_dict=None, format_wildcards=None):
        self.rule = rule
        self.dag = dag

        self.wildcards_dict = wildcards_dict
        self.wildcards = Wildcards(fromdict=self.wildcards_dict)
        self._format_wildcards = (self.wildcards if format_wildcards is None
                                  else Wildcards(fromdict=format_wildcards))

        self.input, input_mapping, self.dependencies = self.rule.expand_input(self.wildcards_dict)
        self.output, output_mapping = self.rule.expand_output(self.wildcards_dict)
        # other properties are lazy to be able to use additional parameters and check already existing files
        self._params = None
        self._log = None
        self._benchmark = None
        self._resources = None
        self._conda_env_file = None
        self._conda_env = None

        self.shadow_dir = None
        self._inputsize = None

        self.restart_times = self.rule.restart_times

        self.dynamic_output, self.dynamic_input = set(), set()
        self.temp_output, self.protected_output = set(), set()
        self.touch_output = set()
        self.subworkflow_input = dict()
        for f in self.output:
            f_ = output_mapping[f]
            if f_ in self.rule.dynamic_output:
                self.dynamic_output.add(f)
            if f_ in self.rule.temp_output:
                self.temp_output.add(f)
            if f_ in self.rule.protected_output:
                self.protected_output.add(f)
            if f_ in self.rule.touch_output:
                self.touch_output.add(f)
        for f in self.input:
            f_ = input_mapping[f]
            if f_ in self.rule.dynamic_input:
                self.dynamic_input.add(f)
            if f_ in self.rule.subworkflow_input:
                self.subworkflow_input[f] = self.rule.subworkflow_input[f_]
            elif "subworkflow" in f.flags:
                self.subworkflow_input[f] = f.flags["subworkflow"]
        self._hash = self.rule.__hash__()
        for o in self.output:
            self._hash ^= o.__hash__()

    def is_valid(self):
        """Check if job is valid"""
        # these properties have to work in dry-run as well. Hence we check them here:
        resources = self.rule.expand_resources(self.wildcards_dict, self.input)
        self.rule.expand_params(self.wildcards_dict, self.input, self.output, resources)
        self.rule.expand_benchmark(self.wildcards_dict)
        self.rule.expand_log(self.wildcards_dict)

    def outputs_older_than_script(self):
        """return output that's older than script, i.e. script has changed"""
        if not self.is_script:
            return
        assert os.path.exists(self.rule.script)# to make sure lstat works
        script_mtime = lstat(self.rule.script).st_mtime
        for f in self.expanded_output:
            if f.exists:
                if not f.is_newer(script_mtime):
                    yield f

    @property
    def threads(self):
        return self.resources._cores

    @property
    def params(self):
        if self._params is None:
            self._params = self.rule.expand_params(self.wildcards_dict,
                                                   self.input,
                                                   self.output,
                                                   self.resources)
        return self._params

    @property
    def log(self):
        if self._log is None:
            self._log = self.rule.expand_log(self.wildcards_dict)
        return self._log

    @property
    def benchmark(self):
        if self._benchmark is None:
            self._benchmark = self.rule.expand_benchmark(self.wildcards_dict)
        return self._benchmark

    @property
    def resources(self):
        if self._resources is None:
            self._resources = self.rule.expand_resources(self.wildcards_dict,
                                                         self.input)
        return self._resources

    @property
    def conda_env_file(self):
        if self._conda_env_file is None:
            expanded_env = self.rule.expand_conda_env(self.wildcards_dict)
            scheme, _, path, *_ = urlparse(expanded_env)
            # Normalize 'file:///my/path.yml' to '/my/path.yml'
            if scheme == 'file' or not scheme:
                self._conda_env_file = path
            else:
                self._conda_env_file = expanded_env
        return self._conda_env_file

    @property
    def conda_env(self):
        if self.conda_env_file:
            if self._conda_env is None:
                self._conda_env = self.dag.conda_envs.get(self.conda_env_file)
            logger.debug("Accessing conda environment {}.".format(self._conda_env))
            if self._conda_env is None:
                raise ValueError("Conda environment {} not found in DAG.".format(self.conda_env_file))
            return self._conda_env.path
        return None

    def archive_conda_env(self):
        """Archive a conda environment into a custom local channel."""
        if self.conda_env_file:
            return self.conda_env.create_archive()
        return None

    @property
    def is_shadow(self):
        return self.rule.shadow_depth is not None

    @property
    def priority(self):
        return self.dag.priority(self)

    @property
    def b64id(self):
        return base64.b64encode((self.rule.name + "".join(self.output)).encode(
            "utf-8")).decode("utf-8")

    @property
    def inputsize(self):
        """
        Return the size of the input files.
        Input files need to be present.
        """
        if self._inputsize is None:
            self._inputsize = sum(f.size for f in self.input)
        return self._inputsize

    @property
    def message(self):
        """ Return the message for this job. """
        try:
            return (self.format_wildcards(self.rule.message) if
                    self.rule.message else None)
        except AttributeError as ex:
            raise RuleException(str(ex), rule=self.rule)
        except KeyError as ex:
            raise RuleException("Unknown variable in message "
                                "of shell command: {}".format(str(ex)),
                                rule=self.rule)

    @property
    def shellcmd(self):
        """ Return the shell command. """
        try:
            return (self.format_wildcards(self.rule.shellcmd) if
                    self.rule.shellcmd else None)
        except AttributeError as ex:
            raise RuleException(str(ex), rule=self.rule)
        except KeyError as ex:
            raise RuleException("Unknown variable when printing "
                                "shell command: {}".format(str(ex)),
                                rule=self.rule)

    @property
    def is_shell(self):
        return self.rule.shellcmd is not None

    @property
    def is_norun(self):
        return self.rule.norun

    @property
    def is_script(self):
        return self.rule.script is not None

    @property
    def is_wrapper(self):
        return self.rule.wrapper is not None

    @property
    def expanded_output(self):
        """ Iterate over output files while dynamic output is expanded. """
        for f, f_ in zip(self.output, self.rule.output):
            if f in self.dynamic_output:
                expansion = self.expand_dynamic(f_)
                if not expansion:
                    yield f_
                for f, _ in expansion:
                    file_to_yield = IOFile(f, self.rule)
                    file_to_yield.clone_flags(f_)
                    yield file_to_yield
            else:
                yield f

    def shadowed_path(self, f):
        """ Get the shadowed path of IOFile f. """
        if not self.shadow_dir:
            return f
        f_ = IOFile(os.path.join(self.shadow_dir, f), self.rule)
        f_.clone_flags(f)
        return f_

    @property
    def dynamic_wildcards(self):
        """ Return all wildcard values determined from dynamic output. """
        combinations = set()
        for f, f_ in zip(self.output, self.rule.output):
            if f in self.dynamic_output:
                for f, w in self.expand_dynamic(f_):
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
        return set(f
                   for f in self.input
                   if not f.exists and not f in self.subworkflow_input)

    @property
    def existing_remote_input(self):
        files = set()

        for f in self.input:
            if f.is_remote:
                if f.exists_remote:
                    files.add(f)
        return files

    @property
    def existing_remote_output(self):
        files = set()

        for f in self.remote_output:
            if f.exists_remote:
                files.add(f)
        return files

    @property
    def missing_remote_input(self):
        return self.remote_input - self.existing_remote_input

    @property
    def missing_remote_output(self):
        return self.remote_output - self.existing_remote_output

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
    def output_mintime_local(self):
        existing = [f.mtime_local for f in self.expanded_output if f.exists]
        if self.benchmark and self.benchmark.exists:
            existing.append(self.benchmark.mtime_local)
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
        if self.benchmark and (requested is None or
                               self.benchmark in requested):
            if not self.benchmark.exists:
                files.add(self.benchmark)

        for f, f_ in zip(self.output, self.rule.output):
            if requested is None or f in requested:
                if f in self.dynamic_output:
                    if not self.expand_dynamic(f_):
                        files.add("{} (dynamic)".format(f_))
                elif not f.exists:
                    files.add(f)
        return files

    @property
    def local_input(self):
        for f in self.input:
            if not f.is_remote:
                yield f

    @property
    def local_output(self):
        for f in self.output:
            if not f.is_remote:
                yield f

    @property
    def remote_input(self):
        for f in self.input:
            if f.is_remote:
                yield f

    @property
    def remote_output(self):
        for f in self.output:
            if f.is_remote:
                yield f

    @property
    def remote_input_newer_than_local(self):
        files = set()
        for f in self.remote_input:
            if (f.exists_remote and f.exists_local) and (
                    f.mtime > f.mtime_local):
                files.add(f)
        return files

    @property
    def remote_input_older_than_local(self):
        files = set()
        for f in self.remote_input:
            if (f.exists_remote and f.exists_local) and (
                    f.mtime < f.mtime_local):
                files.add(f)
        return files

    @property
    def remote_output_newer_than_local(self):
        files = set()
        for f in self.remote_output:
            if (f.exists_remote and f.exists_local) and (
                    f.mtime > f.mtime_local):
                files.add(f)
        return files

    @property
    def remote_output_older_than_local(self):
        files = set()
        for f in self.remote_output:
            if (f.exists_remote and f.exists_local) and (
                    f.mtime < f.mtime_local):
                files.add(f)
        return files

    @property
    def files_to_download(self):
        toDownload = set()

        for f in self.input:
            if f.is_remote:
                if not f.exists_local and f.exists_remote:
                    toDownload.add(f)

        toDownload = toDownload | self.remote_input_newer_than_local
        return toDownload

    @property
    def files_to_upload(self):
        return self.missing_remote_input & self.remote_input_older_than_local

    @property
    def existing_output(self):
        return filter(lambda f: f.exists, self.expanded_output)

    def check_protected_output(self):
        protected = list(filter(lambda f: f.protected, self.expanded_output))
        if protected:
            raise ProtectedOutputException(self.rule, protected)

    def remove_existing_output(self):
        """Clean up both dynamic and regular output before rules actually run
        """
        if self.dynamic_output:
            for f, _ in chain(*map(self.expand_dynamic,
                                   self.rule.dynamic_output)):
                os.remove(f)

        for f, f_ in zip(self.output, self.rule.output):
            try:
                f.remove(remove_non_empty_dir=False)
            except FileNotFoundError:
                #No file == no problem
                pass

    def download_remote_input(self):
        for f in self.files_to_download:
            f.download_from_remote()

    def prepare(self):
        """
        Prepare execution of job.
        This includes creation of directories and deletion of previously
        created dynamic files.
        Creates a shadow directory for the job if specified.
        """

        self.check_protected_output()

        unexpected_output = self.dag.reason(self).missing_output.intersection(
            self.existing_output)
        if unexpected_output:
            logger.warning(
                "Warning: the following output files of rule {} were not "
                "present when the DAG was created:\n{}".format(
                    self.rule, unexpected_output))

        self.remove_existing_output()

        for f, f_ in zip(self.output, self.rule.output):
            f.prepare()

        self.download_remote_input()

        for f in self.log:
            f.prepare()
        if self.benchmark:
            self.benchmark.prepare()

        if not self.is_shadow:
            return
        # Create shadow directory structure
        self.shadow_dir = tempfile.mkdtemp(
            dir=self.rule.workflow.persistence.shadow_path)
        cwd = os.getcwd()
        # Shallow simply symlink everything in the working directory.
        if self.rule.shadow_depth == "shallow":
            for source in os.listdir(cwd):
                link = os.path.join(self.shadow_dir, source)
                os.symlink(os.path.abspath(source), link)
        elif self.rule.shadow_depth == "full":
            snakemake_dir = os.path.join(cwd, ".snakemake")
            for dirpath, dirnames, filenames in os.walk(cwd):
                # Must exclude .snakemake and its children to avoid infinite
                # loop of symlinks.
                if os.path.commonprefix([snakemake_dir, dirpath
                                         ]) == snakemake_dir:
                    continue
                for dirname in dirnames:
                    if dirname == ".snakemake":
                        continue
                    relative_source = os.path.relpath(os.path.join(dirpath,
                                                                   dirname))
                    shadow = os.path.join(self.shadow_dir, relative_source)
                    os.mkdir(shadow)

                for filename in filenames:
                    source = os.path.join(dirpath, filename)
                    relative_source = os.path.relpath(source)
                    link = os.path.join(self.shadow_dir, relative_source)
                    os.symlink(source, link)

    def close_remote(self):
        for f in (self.input + self.output):
            if f.is_remote:
                f.remote_object.close()

    def cleanup(self):
        """ Cleanup output files. """
        to_remove = [f for f in self.expanded_output if f.exists]

        to_remove.extend([f for f in self.remote_input if f.exists_local])
        to_remove.extend([
            f for f in self.remote_output
            if (f.exists_remote if (f.is_remote and f.should_stay_on_remote) else f.exists_local)
        ])
        if to_remove:
            logger.info("Removing output files of failed job {}"
                        " since they might be corrupted:\n{}".format(
                            self, ", ".join(to_remove)))
            for f in to_remove:
                f.remove()

            self.rmdir_empty_remote_dirs()

    @property
    def empty_remote_dirs(self):
        for f in (set(self.output) | set(self.input)):
            if f.is_remote and not f.should_stay_on_remote:
                if os.path.exists(os.path.dirname(f)) and not len(os.listdir(
                        os.path.dirname(f))):
                    yield os.path.dirname(f)

    def rmdir_empty_remote_dirs(self):
        for d in self.empty_remote_dirs:
            try:
                os.removedirs(d)
            except:
                pass  # it's ok if we can't remove the leaf

    def format_wildcards(self, string, **variables):
        """ Format a string with variables from the job. """
        _variables = dict()
        _variables.update(self.rule.workflow.globals)
        _variables.update(dict(input=self.input,
                               output=self.output,
                               params=self.params,
                               wildcards=self._format_wildcards,
                               threads=self.threads,
                               resources=self.resources,
                               log=self.log,
                               version=self.rule.version,
                               rule=self.rule.name, ))
        _variables.update(variables)
        try:
            return format(string, **_variables)
        except NameError as ex:
            raise RuleException("NameError: " + str(ex), rule=self.rule)
        except IndexError as ex:
            raise RuleException("IndexError: " + str(ex), rule=self.rule)

    def properties(self,
                   omit_resources="_cores _nodes".split(),
                   **aux_properties):
        resources = {
            name: res
            for name, res in self.resources.items()
            if name not in omit_resources
        }
        params = {name: value for name, value in self.params.items()}
        properties = {
            "rule": self.rule.name,
            "local": self.dag.workflow.is_local(self.rule),
            "input": self.input,
            "output": self.output,
            "wildcards": self.wildcards,
            "params": params,
            "log": self.log,
            "threads": self.threads,
            "resources": resources,
            "jobid": self.dag.jobid(self)
        }
        properties.update(aux_properties)
        return properties

    def __repr__(self):
        return self.rule.name

    def __eq__(self, other):
        if other is None:
            return False
        return (self.rule == other.rule and
                (self.dynamic_output or
                 self.wildcards_dict == other.wildcards_dict) and
                (self.dynamic_input or self.input == other.input))

    def __lt__(self, other):
        return self.rule.__lt__(other.rule)

    def __gt__(self, other):
        return self.rule.__gt__(other.rule)

    def __hash__(self):
        return self._hash

    def expand_dynamic(self, pattern):
        """ Expand dynamic files. """
        return list(listfiles(pattern,
                              restriction=self.wildcards,
                              omit_value=DYNAMIC_FILL))


class Reason:

    __slots__ = ["_updated_input", "_updated_input_run", "_missing_output",
                 "_incomplete_output", "forced", "noio", "nooutput", "derived"]

    def __init__(self):
        self._updated_input = None
        self._updated_input_run = None
        self._missing_output = None
        self._incomplete_output = None
        self.forced = False
        self.noio = False
        self.nooutput = False
        self.derived = True

    @lazy_property
    def updated_input(self):
        return set()

    @lazy_property
    def updated_input_run(self):
        return set()

    @lazy_property
    def missing_output(self):
        return set()

    @lazy_property
    def incomplete_output(self):
        return set()

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
                if self._missing_output:
                    s.append("Missing output files: {}".format(", ".join(
                        self.missing_output)))
                if self._incomplete_output:
                    s.append("Incomplete output files: {}".format(", ".join(
                        self.incomplete_output)))
                if self._updated_input:
                    updated_input = self.updated_input - self.updated_input_run
                    s.append("Updated input files: {}".format(", ".join(
                        updated_input)))
                if self._updated_input_run:
                    s.append("Input files updated by another job: {}".format(
                        ", ".join(self.updated_input_run)))
        s = "; ".join(s)
        return s

    def __bool__(self):
        return bool(self.updated_input or self.missing_output or self.forced or
                    self.updated_input_run or self.noio or self.nooutput)
