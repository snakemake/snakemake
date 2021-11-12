__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import enum
import os
import sys
import base64
import tempfile
import json
import shutil
import copy

from collections import defaultdict, OrderedDict
from itertools import chain, filterfalse, groupby
from operator import attrgetter, itemgetter

from snakemake.io import (
    IOFile,
    Wildcards,
    Resources,
    _IOFile,
    is_flagged,
    get_flag_value,
    wait_for_files,
)
from snakemake.utils import format, listfiles
from snakemake.exceptions import RuleException, ProtectedOutputException, WorkflowError

from snakemake.logging import logger
from snakemake.common import (
    DYNAMIC_FILL,
    is_local_file,
    parse_uri,
    lazy_property,
    get_uuid,
    TBDString,
)


def format_files(job, io, dynamicio):
    for f in io:
        if f in dynamicio:
            yield f"{f.format_dynamic()} (dynamic)"
        elif is_flagged(f, "pipe"):
            yield f"{f} (pipe)"
        elif is_flagged(f, "service"):
            yield f"{f} (service)"
        elif is_flagged(f, "checkpoint_target"):
            yield TBDString()
        elif is_flagged(f, "sourcecache_entry"):
            orig_path_or_uri = get_flag_value(f, "sourcecache_entry")
            yield f"{orig_path_or_uri} (cached)"
        else:
            yield f


def jobfiles(jobs, type):
    return chain(*map(attrgetter(type), jobs))


class AbstractJob:
    def is_group(self):
        raise NotImplementedError()

    def log_info(self, skip_dynamic=False):
        raise NotImplementedError()

    def log_error(self, msg=None, **kwargs):
        raise NotImplementedError()

    def remove_existing_output(self):
        raise NotImplementedError()

    def download_remote_input(self):
        raise NotImplementedError()

    def properties(self, omit_resources=["_cores", "_nodes"], **aux_properties):
        raise NotImplementedError()

    def reset_params_and_resources(self):
        raise NotImplementedError()


class JobFactory:
    def __init__(self):
        self.cache = dict()

    def new(
        self,
        rule,
        dag,
        wildcards_dict=None,
        format_wildcards=None,
        targetfile=None,
        update=False,
        groupid=None,
    ):
        if rule.is_branched:
            # for distinguishing branched rules, we need input and output in addition
            key = (
                rule.name,
                *rule.output,
                *rule.input,
                *sorted(wildcards_dict.items()),
            )
        else:
            key = (rule.name, *sorted(wildcards_dict.items()))
        if update:
            # cache entry has to be replaced because job shall be constructed from scratch
            obj = Job(rule, dag, wildcards_dict, format_wildcards, targetfile, groupid)
            self.cache[key] = obj
        else:
            try:
                # try to get job from cache
                obj = self.cache[key]
            except KeyError:
                obj = Job(rule, dag, wildcards_dict, format_wildcards, targetfile)
                self.cache[key] = obj
        return obj


class Job(AbstractJob):
    HIGHEST_PRIORITY = sys.maxsize

    obj_cache = dict()

    __slots__ = [
        "rule",
        "dag",
        "wildcards_dict",
        "wildcards",
        "_format_wildcards",
        "input",
        "dependencies",
        "output",
        "_params",
        "_log",
        "_benchmark",
        "_resources",
        "_conda_env_file",
        "_conda_env",
        "shadow_dir",
        "_inputsize",
        "dynamic_output",
        "dynamic_input",
        "temp_output",
        "protected_output",
        "touch_output",
        "subworkflow_input",
        "_hash",
        "_attempt",
        "_group",
        "targetfile",
        "incomplete_input_expand",
        "_params_and_resources_resetted",
    ]

    def __init__(
        self,
        rule,
        dag,
        wildcards_dict=None,
        format_wildcards=None,
        targetfile=None,
        groupid=None,
    ):
        self.rule = rule
        self.dag = dag

        # the targetfile that led to the job
        # it is important to record this, since we need it to submit the
        # job on a cluster. In contrast, an arbitrary targetfile could
        # lead to a different composition of wildcard values (in case of
        # ambiguity in matching).
        self.targetfile = targetfile
        self.wildcards_dict = wildcards_dict
        self.wildcards = Wildcards(fromdict=self.wildcards_dict)
        self._format_wildcards = (
            self.wildcards
            if format_wildcards is None
            else Wildcards(fromdict=format_wildcards)
        )

        (
            self.input,
            input_mapping,
            self.dependencies,
            self.incomplete_input_expand,
        ) = self.rule.expand_input(
            self.wildcards_dict,
            groupid=groupid,
        )

        self.output, output_mapping = self.rule.expand_output(self.wildcards_dict)
        # other properties are lazy to be able to use additional parameters and check already existing files
        self._params = None
        self._log = None
        self._benchmark = None
        self._resources = None
        self._conda_env_spec = None
        self._local_resources = None
        self._conda_env = None
        self._group = None

        # pipe_group will only be set if the job generates or consumes a pipe
        self.pipe_group = None

        self.shadow_dir = None
        self._inputsize = None
        self.is_updated = False
        self._params_and_resources_resetted = False

        self._attempt = self.dag.workflow.attempt

        # TODO get rid of these
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
                sub = f.flags["subworkflow"]
                if f in self.subworkflow_input:
                    other = self.subworkflow_input[f]
                    if sub != other:
                        raise WorkflowError(
                            "The input file {} is ambiguously "
                            "associated with two subworkflows {} "
                            "and {}.".format(f, sub, other),
                            rule=self.rule,
                        )
                self.subworkflow_input[f] = sub

    def updated(self):
        group = self.dag.get_job_group(self)
        groupid = None
        if group is None:
            if self.dag.workflow.run_local or self.is_local:
                groupid = self.dag.workflow.local_groupid
        else:
            groupid = group.jobid

        job = self.dag.job_factory.new(
            self.rule,
            self.dag,
            wildcards_dict=self.wildcards_dict,
            targetfile=self.targetfile,
            update=True,
            groupid=groupid,
        )
        job.is_updated = True
        return job

    def is_valid(self):
        """Check if job is valid"""
        # these properties have to work in dry-run as well. Hence we check them here:
        self.rule.expand_benchmark(self.wildcards_dict)
        self.rule.expand_log(self.wildcards_dict)

    def outputs_older_than_script_or_notebook(self):
        """return output that's older than script, i.e. script has changed"""
        path = self.rule.script or self.rule.notebook
        if not path:
            return
        if self.rule.basedir:
            # needed if rule is included from another subdirectory
            path = self.rule.basedir.join(path).get_path_or_uri()
        if is_local_file(path) and os.path.exists(path):
            script_mtime = os.lstat(path).st_mtime
            for f in self.expanded_output:
                if f.exists:
                    if not f.is_newer(script_mtime):
                        yield f
        # TODO also handle remote file case here.

    @property
    def threads(self):
        return self.resources._cores

    @property
    def params(self):
        if self._params is None:
            self._params = self.rule.expand_params(
                self.wildcards_dict, self.input, self.output, self.resources
            )
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
    def benchmark_repeats(self):
        if self.benchmark is not None:
            return get_flag_value(self.benchmark, "repeat") or 1

    @property
    def group(self):
        if self._group is None:
            self._group = self.rule.expand_group(self.wildcards_dict)
        return self._group

    @group.setter
    def group(self, group):
        self._group = group

    @property
    def attempt(self):
        return self._attempt

    @attempt.setter
    def attempt(self, attempt):
        # reset resources
        self._resources = None
        self._attempt = attempt

    @property
    def resources(self):
        if self._resources is None:
            self._resources = self.rule.expand_resources(
                self.wildcards_dict, self.input, self.attempt
            )
        return self._resources

    @property
    def local_resources(self):
        if self._local_resources is None:
            if self.is_local:
                self._local_resources = self.resources
            else:
                self._local_resources = Resources(
                    fromdict={
                        k: self.resources[k]
                        for k in {"_cores", "_nodes"} & {*self.resources.keys()}
                    }
                )
        return self._local_resources

    def reset_params_and_resources(self):
        if not self._params_and_resources_resetted:
            self._resources = None
            self._params = None
            self._params_and_resources_resetted = True

    @property
    def conda_env_spec(self):
        if self._conda_env_spec is None:
            self._conda_env_spec = self.rule.expand_conda_env(
                self.wildcards_dict, self.params, self.input
            )
        return self._conda_env_spec

    @property
    def conda_env(self):
        if self.conda_env_spec:
            if self._conda_env is None:
                self._conda_env = self.dag.conda_envs.get(
                    (self.conda_env_spec, self.container_img_url)
                )
            return self._conda_env
        return None

    def archive_conda_env(self):
        """Archive a conda environment into a custom local channel."""
        if self.conda_env_spec:
            if self.conda_env.is_named:
                raise WorkflowError(
                    "Workflow archives cannot be created for workflows using named conda environments."
                    "Please use paths to YAML files for all your conda directives.",
                    rule=self.rule,
                )
            return self.conda_env.create_archive()
        return None

    @property
    def needs_singularity(self):
        return self.container_img is not None

    @property
    def container_img_url(self):
        return self.rule.container_img

    @property
    def is_containerized(self):
        return self.rule.is_containerized

    @property
    def container_img(self):
        if self.dag.workflow.use_singularity and self.container_img_url:
            return self.dag.container_imgs[self.container_img_url]
        return None

    @property
    def env_modules(self):
        return self.rule.env_modules

    @property
    def container_img_path(self):
        return self.container_img.path if self.container_img else None

    @property
    def is_shadow(self):
        return self.rule.shadow_depth is not None

    @property
    def priority(self):
        return self.dag.priority(self)

    @property
    def b64id(self):
        return base64.b64encode(
            (self.rule.name + "".join(self.output)).encode("utf-8")
        ).decode("utf-8")

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
        """Return the message for this job."""
        try:
            return (
                self.format_wildcards(self.rule.message) if self.rule.message else None
            )
        except AttributeError as ex:
            raise RuleException(str(ex), rule=self.rule)
        except KeyError as ex:
            raise RuleException(
                "Unknown variable in message " "of shell command: {}".format(str(ex)),
                rule=self.rule,
            )

    @property
    def shellcmd(self):
        """Return the shell command."""
        try:
            return (
                self.format_wildcards(self.rule.shellcmd)
                if self.rule.shellcmd
                else None
            )
        except AttributeError as ex:
            raise RuleException(str(ex), rule=self.rule)
        except KeyError as ex:
            raise RuleException(
                "Unknown variable when printing " "shell command: {}".format(str(ex)),
                rule=self.rule,
            )

    @property
    def is_shell(self):
        return self.rule.is_shell

    @property
    def is_norun(self):
        return self.rule.norun

    @property
    def is_script(self):
        return self.rule.is_script

    @property
    def is_notebook(self):
        return self.rule.is_notebook

    @property
    def is_wrapper(self):
        return self.rule.is_wrapper

    @property
    def is_cwl(self):
        return self.rule.is_cwl

    @property
    def is_template_engine(self):
        return self.rule.is_template_engine

    @property
    def is_run(self):
        return not (
            self.is_shell
            or self.is_norun
            or self.is_script
            or self.is_notebook
            or self.is_wrapper
            or self.is_cwl
        )

    @property
    def is_pipe(self):
        return any(is_flagged(o, "pipe") for o in self.output)

    @property
    def is_service(self):
        return any(is_flagged(o, "service") for o in self.output)

    @property
    def expanded_output(self):
        """Iterate over output files while dynamic output is expanded."""
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
        """Get the shadowed path of IOFile f."""
        if not self.shadow_dir:
            return f
        f_ = IOFile(os.path.join(self.shadow_dir, f), self.rule)
        f_.clone_flags(f)
        return f_

    @property
    def dynamic_wildcards(self):
        """Return all wildcard values determined from dynamic output."""
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
        """Return missing input files."""
        # omit file if it comes from a subworkflow
        return set(
            f for f in self.input if not f.exists and not f in self.subworkflow_input
        )

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
        """Return oldest output file."""
        try:
            mintime = min(
                f.mtime.local_or_remote() for f in self.expanded_output if f.exists
            )
        except ValueError:
            # no existing output
            mintime = None

        if self.benchmark and self.benchmark.exists:
            mintime_benchmark = self.benchmark.mtime.local_or_remote()
            if mintime is not None:
                return min(mintime, mintime_benchmark)
            else:
                return mintime_benchmark

        return mintime

    def missing_output(self, requested):
        def handle_file(f):
            # pipe or service output is always declared as missing
            # (even if it might be present on disk for some reason)
            if is_flagged(f, "pipe") or is_flagged(f, "service") or not f.exists:
                yield f

        if self.dynamic_output:
            for f, f_ in zip(self.output, self.rule.output):
                if f in requested:
                    if f in self.dynamic_output:
                        if not self.expand_dynamic(f_):
                            yield "{} (dynamic)".format(f_)
                    else:
                        yield from handle_file(f)
        else:
            for f in requested:
                yield from handle_file(f)

    @property
    def local_input(self):
        for f in self.input:
            if not f.is_remote:
                yield f

    @property
    def unique_input(self):
        seen = set()

        for element in filterfalse(seen.__contains__, self.input):
            seen.add(element)
            yield element

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
                f.mtime.remote() > f.mtime.local(follow_symlinks=True)
            ):
                files.add(f)
        return files

    @property
    def remote_input_older_than_local(self):
        files = set()
        for f in self.remote_input:
            if (f.exists_remote and f.exists_local) and (
                f.mtime.remote() < f.mtime.local(follow_symlinks=True)
            ):
                files.add(f)
        return files

    @property
    def remote_output_newer_than_local(self):
        files = set()
        for f in self.remote_output:
            if (f.exists_remote and f.exists_local) and (
                f.mtime.remote() > f.mtime.local(follow_symlinks=True)
            ):
                files.add(f)
        return files

    @property
    def remote_output_older_than_local(self):
        files = set()
        for f in self.remote_output:
            if (f.exists_remote and f.exists_local) and (
                f.mtime.remote() < f.mtime.local(follow_symlinks=True)
            ):
                files.add(f)
        return files

    @property
    def files_to_download(self):
        toDownload = set()

        for f in self.input:
            if f.is_remote:
                if (not f.exists_local and f.exists_remote) and (
                    not self.rule.norun or f.remote_object.keep_local
                ):
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
            raise ProtectedOutputException(self, protected)

    def remove_existing_output(self):
        """Clean up both dynamic and regular output before rules actually run"""
        if self.dynamic_output:
            for f, _ in chain(*map(self.expand_dynamic, self.rule.dynamic_output)):
                os.remove(f)

        for f, f_ in zip(self.output, self.rule.output):
            try:
                # remove_non_empty_dir only applies to directories which aren't
                # flagged with directory().
                f.remove(remove_non_empty_dir=False)
            except FileNotFoundError:
                # No file == no problem
                pass

        for f in self.log:
            f.remove(remove_non_empty_dir=False)

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
            self.existing_output
        )
        if unexpected_output:
            logger.warning(
                "Warning: the following output files of rule {} were not "
                "present when the DAG was created:\n{}".format(
                    self.rule, unexpected_output
                )
            )

        self.remove_existing_output()

        # Create tmpdir if necessary
        if self.resources.get("tmpdir"):
            os.makedirs(self.resources.tmpdir, exist_ok=True)

        for f, f_ in zip(self.output, self.rule.output):
            f.prepare()

        self.download_remote_input()

        for f in self.log:
            f.prepare()
        if self.benchmark:
            self.benchmark.prepare()

        # wait for input files
        wait_for_files(self.input, latency_wait=self.dag.workflow.latency_wait)

        if not self.is_shadow:
            return

        # Create shadow directory structure
        self.shadow_dir = tempfile.mkdtemp(
            dir=self.rule.workflow.persistence.shadow_path
        )
        cwd = os.getcwd()

        # "minimal" creates symlinks only to the input files in the shadow directory
        # "copy-minimal" creates copies instead
        if (
            self.rule.shadow_depth == "minimal"
            or self.rule.shadow_depth == "copy-minimal"
        ):
            # Re-create the directory structure in the shadow directory
            for (f, d) in set(
                [
                    (item, os.path.dirname(item))
                    for sublist in [self.input, self.output, self.log]
                    if sublist is not None
                    for item in sublist
                ]
            ):
                if d and not os.path.isabs(d):
                    rel_path = os.path.relpath(d)
                    # Only create subdirectories
                    if not rel_path.split(os.path.sep)[0] == "..":
                        os.makedirs(
                            os.path.join(self.shadow_dir, rel_path), exist_ok=True
                        )
                    else:
                        raise RuleException(
                            "The following file name references a parent directory relative to your workdir.\n"
                            'This isn\'t supported for shadow: "{}". Consider using an absolute path instead.\n{}'.format(
                                f, self.rule.shadow_depth
                            ),
                            rule=self.rule,
                        )

            # Symlink or copy the input files
            if self.rule.shadow_depth == "copy-minimal":
                for rel_path in set(
                    [os.path.relpath(f) for f in self.input if not os.path.isabs(f)]
                ):
                    copy = os.path.join(self.shadow_dir, rel_path)
                    shutil.copy(rel_path, copy)
            else:
                for rel_path in set(
                    [os.path.relpath(f) for f in self.input if not os.path.isabs(f)]
                ):
                    link = os.path.join(self.shadow_dir, rel_path)
                    original = os.path.relpath(rel_path, os.path.dirname(link))
                    os.symlink(original, link)

        # Shallow simply symlink everything in the working directory.
        elif self.rule.shadow_depth == "shallow":
            for source in os.listdir(cwd):
                link = os.path.join(self.shadow_dir, source)
                os.symlink(os.path.abspath(source), link)
        elif self.rule.shadow_depth == "full":
            snakemake_dir = os.path.join(cwd, ".snakemake")
            for dirpath, dirnames, filenames in os.walk(cwd):
                # Must exclude .snakemake and its children to avoid infinite
                # loop of symlinks.
                if os.path.commonprefix([snakemake_dir, dirpath]) == snakemake_dir:
                    continue
                for dirname in dirnames:
                    if dirname == ".snakemake":
                        continue
                    relative_source = os.path.relpath(os.path.join(dirpath, dirname))
                    shadow = os.path.join(self.shadow_dir, relative_source)
                    os.mkdir(shadow)

                for filename in filenames:
                    source = os.path.join(dirpath, filename)
                    relative_source = os.path.relpath(source)
                    link = os.path.join(self.shadow_dir, relative_source)
                    os.symlink(source, link)

    def close_remote(self):
        for f in self.input + self.output:
            if f.is_remote:
                f.remote_object.close()

    def cleanup(self):
        """Cleanup output files."""
        to_remove = [f for f in self.expanded_output if f.exists]

        to_remove.extend(
            [
                f
                for f in self.remote_output
                if (
                    f.exists_remote
                    if (f.is_remote and f.should_stay_on_remote)
                    else f.exists_local
                )
            ]
        )
        if to_remove:
            logger.info(
                "Removing output files of failed job {}"
                " since they might be corrupted:\n{}".format(self, ", ".join(to_remove))
            )
            for f in to_remove:
                f.remove()

    def format_wildcards(self, string, **variables):
        """Format a string with variables from the job."""
        _variables = dict()
        _variables.update(self.rule.workflow.globals)
        _variables.update(
            dict(
                input=self.input,
                output=self.output,
                params=self.params,
                wildcards=self._format_wildcards,
                threads=self.threads,
                resources=self.resources,
                log=self.log,
                jobid=self.jobid,
                version=self.rule.version,
                name=self.name,
                rule=self.rule.name,
                rulename=self.rule.name,
                bench_iteration=None,
            )
        )
        _variables.update(variables)
        try:
            return format(string, **_variables)
        except NameError as ex:
            raise RuleException("NameError: " + str(ex), rule=self.rule)
        except IndexError as ex:
            raise RuleException("IndexError: " + str(ex), rule=self.rule)
        except Exception as ex:
            raise WorkflowError(
                f"Error when formatting '{string}' for rule {self.rule.name}. {ex}"
            )

    def properties(self, omit_resources=["_cores", "_nodes"], **aux_properties):
        resources = {
            name: res
            for name, res in self.resources.items()
            if name not in omit_resources
        }
        params = {name: value for name, value in self.params.items()}
        properties = {
            "type": "single",
            "rule": self.rule.name,
            "local": self.is_local,
            "input": self.input,
            "output": self.output,
            "wildcards": self.wildcards_dict,
            "params": params,
            "log": self.log,
            "threads": self.threads,
            "resources": resources,
            "jobid": self.dag.jobid(self),
        }
        properties.update(aux_properties)

        try:
            return json.dumps(properties)
        except TypeError:
            del properties["params"]
            return json.dumps(properties)

    @property
    def is_local(self):
        return self.dag.workflow.is_local(self.rule)

    def __repr__(self):
        return self.rule.name

    def __lt__(self, other):
        return self.rule.__lt__(other.rule)

    def __gt__(self, other):
        return self.rule.__gt__(other.rule)

    def expand_dynamic(self, pattern):
        """Expand dynamic files."""
        return list(
            listfiles(pattern, restriction=self.wildcards, omit_value=DYNAMIC_FILL)
        )

    def is_group(self):
        return False

    def log_info(self, skip_dynamic=False, indent=False, printshellcmd=True):
        # skip dynamic jobs that will be "executed" only in dryrun mode
        if skip_dynamic and self.dag.dynamic(self):
            return

        priority = self.priority
        logger.job_info(
            jobid=self.dag.jobid(self),
            msg=self.message,
            name=self.rule.name,
            local=self.dag.workflow.is_local(self.rule),
            input=list(format_files(self, self.input, self.dynamic_input)),
            output=list(format_files(self, self.output, self.dynamic_output)),
            log=list(self.log),
            benchmark=self.benchmark,
            wildcards=self.wildcards_dict,
            reason=str(self.dag.reason(self)),
            resources=self.resources,
            priority="highest" if priority == Job.HIGHEST_PRIORITY else priority,
            threads=self.threads,
            indent=indent,
            is_checkpoint=self.rule.is_checkpoint,
            printshellcmd=printshellcmd,
            is_handover=self.rule.is_handover,
        )
        logger.shellcmd(self.shellcmd, indent=indent)

        if self.dynamic_output:
            logger.info(
                "Subsequent jobs will be added dynamically "
                "depending on the output of this job",
                indent=True,
            )

    def log_error(self, msg=None, indent=False, **kwargs):
        logger.job_error(
            name=self.rule.name,
            jobid=self.dag.jobid(self),
            output=list(format_files(self, self.output, self.dynamic_output)),
            log=list(self.log),
            conda_env=self.conda_env.address if self.conda_env else None,
            aux=kwargs,
            indent=indent,
            shellcmd=self.shellcmd,
        )
        if msg is not None:
            logger.error(msg)

    def register(self):
        self.dag.workflow.persistence.started(self)

    def get_wait_for_files(self):
        wait_for_files = []
        wait_for_files.extend(self.local_input)
        wait_for_files.extend(
            f for f in self.remote_input if not f.should_stay_on_remote
        )

        if self.shadow_dir:
            wait_for_files.append(self.shadow_dir)
        if (
            self.dag.workflow.use_conda
            and self.conda_env
            and not self.conda_env.is_named
            and not self.conda_env.is_containerized
        ):
            # Named or containerized envs are not present on the host FS,
            # hence we don't need to wait for them.
            wait_for_files.append(self.conda_env.address)
        return wait_for_files

    @property
    def jobid(self):
        return self.dag.jobid(self)

    def uuid(self):
        return str(
            get_uuid(
                f"{self.rule.name}:{','.join(sorted(f'{w}:{v}' for w, v in self.wildcards_dict.items()))}"
            )
        )

    def postprocess(
        self,
        upload_remote=True,
        handle_log=True,
        handle_touch=True,
        error=False,
        ignore_missing_output=False,
        assume_shared_fs=True,
        latency_wait=None,
        keep_metadata=True,
    ):
        if self.dag.is_edit_notebook_job(self):
            # No postprocessing necessary, we have just created the skeleton notebook and
            # execution will anyway stop afterwards.
            return
        if assume_shared_fs:
            if not error and handle_touch:
                self.dag.handle_touch(self)
            if handle_log:
                self.dag.handle_log(self)
            if not error:
                self.dag.check_and_touch_output(
                    self, wait=latency_wait, ignore_missing_output=ignore_missing_output
                )
            self.dag.unshadow_output(self, only_log=error)
            if not error:
                self.dag.handle_remote(self, upload=upload_remote)
                self.dag.handle_protected(self)
            self.close_remote()
        else:
            if not error:
                self.dag.check_and_touch_output(
                    self, wait=latency_wait, no_touch=True, force_stay_on_remote=True
                )
        if not error:
            try:
                self.dag.workflow.persistence.finished(
                    self, keep_metadata=keep_metadata
                )
            except IOError as e:
                raise WorkflowError(
                    "Error recording metadata for finished job "
                    "({}). Please ensure write permissions for the "
                    "directory {}".format(e, self.dag.workflow.persistence.path)
                )

    @property
    def name(self):
        return self.rule.name

    @property
    def priority(self):
        return self.dag.priority(self)

    @property
    def products(self):
        products = list(self.output)
        if self.benchmark:
            products.append(self.benchmark)
        products.extend(self.log)
        return products

    def get_targets(self):
        return [self.targetfile or self.rule.name]

    @property
    def is_branched(self):
        return self.rule.is_branched

    @property
    def rules(self):
        return [self.rule.name]

    @property
    def restart_times(self):
        return self.rule.restart_times

    @property
    def is_checkpoint(self):
        return self.rule.is_checkpoint

    def __len__(self):
        return 1


class GroupJobFactory:
    def __init__(self):
        self.cache = dict()

    def new(self, id, jobs, resources):
        jobs = frozenset(jobs)
        key = (id, jobs)
        try:
            obj = self.cache[key]
        except KeyError:
            obj = GroupJob(id, jobs, resources)
            self.cache[key] = obj
        return obj


class GroupJob(AbstractJob):

    obj_cache = dict()

    __slots__ = [
        "groupid",
        "jobs",
        "_resources",
        "_input",
        "_output",
        "_log",
        "_inputsize",
        "_all_products",
        "_attempt",
        "toposorted",
        "_jobid",
    ]

    def __init__(self, id, jobs, global_resources):
        self.groupid = id
        self.jobs = jobs
        self.global_resources = global_resources
        self.toposorted = None
        self._resources = None
        self._local_resources = None
        self._input = None
        self._output = None
        self._log = None
        self._inputsize = None
        self._all_products = None
        self._attempt = self.dag.workflow.attempt
        self._jobid = None

    @property
    def dag(self):
        return next(iter(self.jobs)).dag

    def merge(self, other):
        assert other.groupid == self.groupid
        self.jobs = self.jobs | other.jobs

    def finalize(self):
        if self.toposorted is None:
            # Jobs connected by a pipe must be run simultaneously. Thus, we need to
            # trick toposort into putting them on the same level. This is done by
            # causing all jobs in the pipe group to use each other's dependencies.

            # First, we organize all the jobs in the group into a dict according to
            # their pipe_group
            pipe_groups = defaultdict(list)
            for name, group in groupby(self.jobs, attrgetter("pipe_group")):
                if name is not None:
                    pipe_groups[name].extend(group)

            # Then, for each pipe_group, we find the dependencies of every job in the
            # group, filtering out any dependencies that are, themselves, in the group
            pipe_dependencies = {}
            for name, group in pipe_groups.items():
                pipe_dependencies[name] = set(
                    d
                    for job in group
                    for d in self.dag.dependencies[job]
                    if d not in group
                )

            # Collect every job's dependencies into a definitive mapping. Regular jobs
            # get their dependencies directly from dag.dependencies, while pipe jobs
            # get their dependencies from the pipe_dependencies calculated above.
            dependencies = {}
            for job in self.jobs:
                if job.pipe_group in pipe_dependencies:
                    dependencies[job] = pipe_dependencies[job.pipe_group]
                else:
                    dependencies[job] = self.dag.dependencies[job]

            self.toposorted = [*self.toposort(self.jobs, dependencies)]

    def toposort(self, group, dependencies):
        """Sort group of jobs into layers of dependency.

        Args:
            group (Iterable[jobs]): List of jobs to sort
            dependencies (Dict[job, Iterable[jobs]]): Dictionary of job keys mapped to
                a list of all jobs they depend on. Dependencies should have one key
                entry for every job in group.
        """
        from toposort import toposort as _toposort

        # _toposort takes a dict: {job: {dependencies...}, ...}. Here, just filter
        # out all the dependencies that aren't in the group.
        yield from _toposort(
            {job: {dep for dep in dependencies[job] if dep in group} for job in group}
        )

    def __iter__(self):
        if self.toposorted is None:
            yield from self.jobs
            return

        # The scheduler currently relies on Pipe jobs being returned in dependency order
        # (i.e. if a -> b, yield a before b). finalize(), which produces self.toposort,
        # puts all pipe jobs at the same toposort level for resource calculation
        # purposes, so if we yield directly from toposort, we can't guarentee the order.
        # Here, we filter through jobs to find pipe groups, then rerun toposort on each
        # pipe group.
        pipe_groups = defaultdict(set)
        for sibling in self.toposorted:
            for job in sibling:
                if job.pipe_group is None:
                    yield job
                    continue
                pipe_groups[job.pipe_group].add(job)

            for group in pipe_groups.values():
                yield from chain(*self.toposort(group, self.dag.dependencies))

    def __repr__(self):
        return "JobGroup({},{})".format(self.groupid, repr(self.jobs))

    def __contains__(self, job):
        return job in self.jobs

    def is_group(self):
        return True

    @property
    def all_products(self):
        if self._all_products is None:
            self._all_products = set(f for job in self.jobs for f in job.products)
        return self._all_products

    @property
    def is_checkpoint(self):
        return any(job.is_checkpoint for job in self.jobs)

    @property
    def is_updated(self):
        return any(job.is_updated for job in self.jobs)

    def log_info(self, skip_dynamic=False):
        logger.group_info(groupid=self.groupid)
        for job in sorted(self.jobs, key=lambda j: j.rule.name):
            job.log_info(skip_dynamic, indent=True)

    def log_error(self, msg=None, **kwargs):
        logger.group_error(groupid=self.groupid)
        for job in self.jobs:
            job.log_error(msg=msg, indent=True, **kwargs)

    def register(self):
        for job in self.jobs:
            job.register()

    def remove_existing_output(self):
        for job in self.jobs:
            job.remove_existing_output()

    def reset_params_and_resources(self):
        for job in self.jobs:
            job.reset_params_and_resources()

    def download_remote_input(self):
        for job in self.jobs:
            job.download_remote_input()

    def get_wait_for_files(self):
        local_input = [
            f
            for job in self.jobs
            for f in job.local_input
            if f not in self.all_products
        ]
        remote_input = [
            f
            for job in self.jobs
            for f in job.remote_input
            if f not in self.all_products
        ]

        wait_for_files = []
        wait_for_files.extend(local_input)
        wait_for_files.extend(f for f in remote_input if not f.should_stay_on_remote)

        for job in self.jobs:
            if job.shadow_dir:
                wait_for_files.append(job.shadow_dir)
            if (
                self.dag.workflow.use_conda
                and job.conda_env
                and not job.conda_env.is_named
            ):
                wait_for_files.append(job.conda_env.address)
        return wait_for_files

    @property
    def resources(self):
        if self._resources is None:
            self._resources = self._calculate_resources()
        return Resources(fromdict=self._resources)

    @property
    def local_resources(self):
        if self._local_resources is None:
            if self.is_local:
                self._local_resources = self.resources
            else:
                self._local_resources = Resources(
                    fromdict={
                        k: self.resources[k]
                        for k in {"_cores", "_nodes"} & {*self.resources.keys()}
                    }
                )
        return self._local_resources

    @property
    def input(self):
        if self._input is None:
            self._input = [
                f for job in self.jobs for f in job.input if f not in self.all_products
            ]
        return self._input

    @property
    def output(self):
        all_input = set(f for job in self.jobs for f in job.input)
        if self._output is None:
            self._output = [
                f for job in self.jobs for f in job.output if f not in all_input
            ]
        return self._output

    @property
    def log(self):
        if self._log is None:
            self._log = [f for job in self.jobs for f in job.log]
        return self._log

    @property
    def products(self):
        all_input = set(f for job in self.jobs for f in job.input)
        return [f for job in self.jobs for f in job.products if f not in all_input]

    def properties(self, omit_resources=["_cores", "_nodes"], **aux_properties):
        resources = {
            name: res
            for name, res in self.resources.items()
            if name not in omit_resources
        }
        properties = {
            "type": "group",
            "groupid": self.groupid,
            "local": self.is_local,
            "input": self.input,
            "output": self.output,
            "threads": self.threads,
            "resources": resources,
            "jobid": self.jobid,
        }
        properties.update(aux_properties)

        return json.dumps(properties)

    @property
    def jobid(self):
        if not self._jobid:
            # The uuid of the last job is sufficient to uniquely identify the group job.
            # This is true because each job can only occur in one group job.
            # Additionally, this is the most stable id we can get, even if the group
            # changes by adding more upstream jobs, e.g. due to groupid usage in input
            # functions (see Dag.update_incomplete_input_expand_jobs())
            last_job = sorted(self.toposorted[-1])[-1]
            self._jobid = last_job.uuid()
        return self._jobid

    def cleanup(self):
        for job in self.jobs:
            job.cleanup()

    def postprocess(self, error=False, **kwargs):
        for job in self.jobs:
            job.postprocess(error=error, **kwargs)
        # remove all pipe and service outputs since all jobs of this group are done and the
        # outputs are no longer needed
        for job in self.jobs:
            for f in job.output:
                if is_flagged(f, "pipe") or is_flagged(f, "service"):
                    f.remove()

    @property
    def name(self):
        return str(self.groupid)

    def check_protected_output(self):
        for job in self.jobs:
            job.check_protected_output()

    @property
    def dynamic_input(self):
        return [
            f
            for job in self.jobs
            for f in job.dynamic_input
            if f not in self.all_products
        ]

    @property
    def inputsize(self):
        if self._inputsize is None:
            self._inputsize = sum(f.size for f in self.input)
        return self._inputsize

    @property
    def priority(self):
        return max(self.dag.priority(job) for job in self.jobs)

    @property
    def is_local(self):
        return all(job.is_local for job in self.jobs)

    def format_wildcards(self, string, **variables):
        """Format a string with variables from the job."""
        _variables = dict()
        _variables.update(self.dag.workflow.globals)
        _variables.update(
            dict(
                input=self.input,
                output=self.output,
                threads=self.threads,
                jobid=self.jobid,
                name=self.name,
                rule="GROUP",
                rulename="GROUP",
                resources=self.resources,
            )
        )
        _variables.update(variables)
        try:
            return format(string, **_variables)
        except NameError as ex:
            raise WorkflowError(
                "NameError with group job {}: {}".format(self.jobid, str(ex))
            )
        except IndexError as ex:
            raise WorkflowError(
                "IndexError with group job {}: {}".format(self.jobid, str(ex))
            )
        except Exception as ex:
            raise WorkflowError(
                f"Error when formatting {string} for group job {self.jobid}: {ex}"
            )

    @property
    def threads(self):
        return self.resources["_cores"]

    def get_targets(self):
        # jobs without output are targeted by rule name
        targets = [job.rule.name for job in self.jobs if not job.products]
        targets.extend(self.products)
        return targets

    @property
    def attempt(self):
        return self._attempt

    @attempt.setter
    def attempt(self, attempt):
        # reset resources
        self._resources = None
        self._attempt = attempt

    @property
    def is_branched(self):
        return any(job.is_branched for job in self.jobs)

    @property
    def needs_singularity(self):
        return any(job.needs_singularity for job in self.jobs)

    @property
    def rules(self):
        return [job.rule.name for job in self.jobs]

    @property
    def expanded_output(self):
        """Yields the entire expanded output of all jobs"""
        for job in self.jobs:
            yield from job.expanded_output

    @property
    def restart_times(self):
        return max(job.restart_times for job in self.jobs)

    def __len__(self):
        return len(self.jobs)

    def __hash__(self):
        return hash(self.jobs)

    def __eq__(self, other):
        if not isinstance(other, AbstractJob):
            return False
        if other.is_group():
            return self.jobs == other.jobs
        else:
            return False

    def _calculate_resources(self):
        total_resources = defaultdict(int)
        total_resources["_nodes"] = 1
        blocks = []
        # iterate over siblings that can be executed in parallel
        for siblings in self.toposorted:
            # Total resource requirements for this toposort layer
            block_resources = {}

            job_resources = []
            pipe_resources = defaultdict(list)
            for job in siblings:
                # Get resources, filtering out FileNotFoundErrors. List items will
                # be job resources objects with (resource: value)
                # [
                #   { "runtime": 5, "threads": 2, "tmpdir": "/tmp" },
                #   { "runtime": 15, "tmpdir": "/tmp"},
                #   ...
                # ]
                # Pipe jobs and regular jobs are put in seperate lists.
                try:
                    # Remove any TBDStrings from values. These will typically arise
                    # here because the default mem_mb and disk_mb are based off of
                    # input file size, and intermediate files in the group are not yet
                    # generated. Thus rules consuming such files must thus explicitely
                    # specify their resources
                    res = {
                        k: res
                        for k, res, in job.resources.items()
                        if not isinstance(res, TBDString)
                    }
                    if job.pipe_group:
                        pipe_resources[job.pipe_group].append(res)
                    else:
                        job_resources.append(res)
                except FileNotFoundError:
                    # Skip job if resource evaluation leads to a file not found error.
                    # This will be caused by an inner job, which needs files created by the same group.
                    # All we can do is to ignore such jobs for now.
                    continue

            # Jobs in pipe groups must be run simultaneously, so we merge all the
            # resources of each pipe group into one big "job". Most resources are
            # summed, except for cores, which are combined via max().
            for pipe in pipe_resources.values():
                job_resources.append(
                    self._simple_merge(pipe, use_max=False, merge_via_other=["_cores"])
                )

            # Set of resource types requested in at least one job
            resource_types = list(set(chain(*job_resources)))
            int_resources = OrderedDict()
            # Sort all integer resources in job_resources into int_resources. Resources
            # defined as a string are placed immediately into block_resources.
            for res in resource_types:
                if res == "_nodes":
                    continue

                # We use 0 as the default. If we ever want to allow None as an allowable
                # "2nd" option for string resources (e.g. tmp path is set on some jobs
                # but not others), this will cause a problem
                values = [resources.get(res, 0) for resources in job_resources]

                if self._is_string_resource(res, values):
                    block_resources[res] = values[0]
                else:
                    int_resources[res] = values

            # Collect values from global_resources to use as constraints.
            constraints = OrderedDict(
                [
                    (name, self.global_resources.get(name, None))
                    for name in int_resources
                ]
            )

            # For now, we are unable to handle a constraint on runtime, so ignore.
            # Jobs requesting too much runtime will still get flagged by the
            # scheduler
            if "runtime" in constraints:
                constraints["runtime"] = None

            layers = self._get_layers(int_resources, constraints.values(), "runtime")

            # We extract runtime out of each layer by finding its position in the job
            # tuple, which will be at the same index as in the OrderedDict
            # constrained_resources.
            # Add the max runtime in each layer across all layers
            if "runtime" in int_resources:
                runtime_i = list(int_resources).index("runtime")

                runtimes = [max(list(zip(*layer))[runtime_i]) for layer in layers]
                block_resources["runtime"] = sum(runtimes)

            # Determine the total amount of each resource by summing up each layer,
            # then taking the max:

            # In each layer, sum across all resource types within the layer,
            # similar to summing along axis 0 in numpy:
            # [
            #   ( 3 ^ , 4 ^ , 1 ),
            #   ( 2 | , 1 | , 6 ),
            #   ( 1 | , 4 | , 0 ),
            # ]
            sums = [[sum(r) for r in zip(*layer)] for layer in layers]

            # Match up layered resources with definitions in int_resources. Runtime has
            # already been calculated above via a different method, so we skip it.
            for i, name in enumerate(int_resources):
                if name != "runtime":
                    block_resources[name] = max([layer[i] for layer in sums])

            blocks.append(block_resources)

        if self.dag.workflow.run_local:
            return Resources(
                fromdict={**self._simple_merge(blocks, use_max=False), "_nodes": 1}
            )
        else:
            return Resources(
                fromdict={
                    **self._simple_merge(
                        blocks, use_max=True, merge_via_other=["runtime"]
                    ),
                    "_nodes": 1,
                }
            )

    def _is_string_resource(self, name, values):
        # If any one of the values provided for a resource is not an int, we
        # can't process it in any way. So we constrain all such resource to be
        # the same
        if all([isinstance(val, int) for val in values]):
            return False
        else:
            unique = set(values)
            if len(unique) > 1:
                raise WorkflowError(
                    'Failed to group jobs in group "{group}". Resource {name} '
                    "is a string but not all group jobs require the same value. "
                    "Observed values: {values}.".format(
                        group=self.groupid, name=name, values=unique
                    )
                )
            return True

    def _simple_merge(self, jobs, skip=[], use_max=True, merge_via_other=[]):
        grouped = {}
        for job in jobs:
            # Wrap every value in job with a list so that lists can be merged later
            job_l = {k: [v] for k, v in job.items()}

            # Merge two dicts together, merging key-values found in both into a
            # list. Code adapted from
            # https://stackoverflow.com/a/11012181/16980632
            grouped = {
                **grouped,
                **job_l,
                **{k: grouped[k] + job_l[k] for k in grouped.keys() & job_l},
            }

        ret = {}
        for res, values in grouped.items():
            if res in skip:
                continue

            merge_via_max = any(
                [
                    use_max and res not in merge_via_other,
                    not use_max and res in merge_via_other,
                ]
            )
            if self._is_string_resource(res, values):
                ret[res] = values[0]
            elif merge_via_max:
                ret[res] = max(values)
            else:
                ret[res] = sum(values)
        return ret

    def _check_constraint(self, resources, constraints):
        sums = [sum(res) for res in zip(*resources)]
        for s, constraint in zip(sums, constraints):
            if constraint:
                layers, mod = divmod(s, constraint)
            else:
                layers = 1
                mod = 0

            # If mod not 0, we add 1 to the number of layers. We then subtract
            # 1, so that if everything fits within the constraint we have 0,
            # otherwise, some number higher than 0. Finally, we convert to bool.
            # If the result is 0 or negative, it fits. If greater, it doesn't
            # fit so we return False
            if bool(max(0, layers + int(bool(mod)) - 1)):
                return False
        return True

    def _get_layers(self, resources, constraints, sortby=None):
        """Calculate required consecutive job layers.

        Layers are used to keep resource requirements within given
        constraint. For instance, if the jobs together require 50 threads,
        but only 12 are available, we will use 5 layers. If multiple constraints are
        used, all will be considered and met. Any constraints given as None will be
        treated as infinite.
        """

        # Calculates the ratio of resource to constraint. E.g, if the resource is 12
        #  cores, and the constraint is 16, it will return 0.75. This is done for
        # every resource type in the group, returning the result in a list
        def _proportion(group):
            return [r / c if c else 0 for r, c in zip(group, constraints)]

        # Return the highest _proportion item in the list
        def _highest_proportion(group):
            return max(_proportion(group))

        rows = [[]]

        # Resources is an OrderedDict: { resource: [val1, val2, val3], ...}.
        # By zipping, we combine the vals into tuples based on job, 1 tuple per
        # job: [ (val1, 1_val1, 2_val1), ...]. In each tuple, the resources
        # will remain in the same order as the OrderedDict, so their identity
        # can be extracted later.
        resource_groups = zip(*resources.values())

        # Sort by _proportion highest to lowest
        pre_sorted = sorted(resource_groups, key=_highest_proportion, reverse=True)

        # If a key is provided (e.g. runtime), we sort again by that key from
        # highest to lowest
        if sortby and sortby in resources:
            # Find the position of the key in the job tuple
            i = list(resources).index(sortby)
            pre_sorted = sorted(pre_sorted, key=itemgetter(i), reverse=True)

        for group in pre_sorted:
            appended = False

            # Check each row for space, starting with the first.
            for row in rows:
                if not appended and self._check_constraint(row + [group], constraints):
                    row.append(group)
                    appended = True

            # If the final "row" in rows has something, we add a new empty
            # row. That way, we guarantee we have a row with space
            if len(rows[-1]) > 0:
                rows.append([])

            # If not appended, that means a rule required more resource
            # than allowed by the constraint. This should only be possible for pipe
            # jobs, which must be run simultaneously.
            if not appended:
                too_high = []
                for i, val in enumerate(_proportion(group)):
                    if val > 1:
                        too_high.append(
                            (list(resources)[i], group[i], list(constraints)[i])
                        )

                error_text = [
                    f"\t{res}: {amount}/{constraint}"
                    for res, amount, constraint in too_high
                ]
                raise WorkflowError(
                    "Not enough resources were provided. This error is typically\n"
                    "caused by a Pipe group requiring too many resources. Note\n"
                    "that all resources including 'runtime' are summed across\n"
                    "every member of the pipe group, except for cores, which is\n"
                    "calculated via max().\n\n"
                    "Excess Resources:\n" + "\n".join(error_text)
                )

        # Remove final empty row. (The above loop ends each cycle by ensuring
        # there's an empty row)
        rows.pop()
        return rows


class Reason:

    __slots__ = [
        "_updated_input",
        "_updated_input_run",
        "_missing_output",
        "_incomplete_output",
        "input_changed",
        "code_changed",
        "params_changed",
        "software_stack_changed",
        "forced",
        "noio",
        "nooutput",
        "derived",
        "pipe",
        "service",
        "target",
        "finished",
        "cleanup_metadata_instructions",
    ]

    def __init__(self):
        self.finished = False
        self._updated_input = None
        self._updated_input_run = None
        self._missing_output = None
        self._incomplete_output = None
        self.params_changed = False
        self.code_changed = False
        self.software_stack_changed = False
        self.input_changed = False
        self.forced = False
        self.noio = False
        self.nooutput = False
        self.derived = True
        self.pipe = False
        self.service = False
        self.cleanup_metadata_instructions = None

    def set_cleanup_metadata_instructions(self, job):
        self.cleanup_metadata_instructions = (
            "To ignore these changes, run snakemake "
            f"--cleanup-metadata {' '.join(job.expanded_output)}"
        )

    def is_provenance_triggered(self):
        """Return True if reason is triggered by provenance information."""
        return (
            self.params_changed
            or self.code_changed
            or self.software_stack_changed
            or self.input_changed
        )

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

    def mark_finished(self):
        "called if the job has been run"
        self.finished = True

    def __str__(self):
        def format_file(f):
            if is_flagged(f, "sourcecache_entry"):
                return f"{get_flag_value(f, 'sourcecache_entry')} (cached)"
            else:
                return f

        def format_files(files):
            return ", ".join(map(format_file, files))

        s = list()
        if self.forced:
            s.append("Forced execution")
        else:
            if self.noio:
                s.append(
                    "Rules with neither input nor " "output files are always executed."
                )
            elif self.nooutput:
                s.append(
                    "Rules with a run or shell declaration but no output "
                    "are always executed."
                )
            else:
                if self._missing_output:
                    s.append(
                        f"Missing output files: {format_files(self.missing_output)}"
                    )
                if self._incomplete_output:
                    s.append(
                        f"Incomplete output files: {format_files(self.incomplete_output)}"
                    )
                if self._updated_input:
                    updated_input = self.updated_input - self.updated_input_run
                    s.append(f"Updated input files: {format_files(updated_input)}")
                if self._updated_input_run:
                    s.append(
                        f"Input files updated by another job: {format_files(self.updated_input_run)}"
                    )
                if self.pipe:
                    s.append(
                        "Output file is a pipe and has to be filled for consuming job."
                    )
                if self.service:
                    s.append(
                        "Job provides a service which has to be kept active until all consumers are finished."
                    )

                if self.input_changed:
                    s.append("Set of input files has changed since last execution")
                if self.code_changed:
                    s.append("Code has changed since last execution")
                if self.params_changed:
                    s.append("Params have changed since last execution")
                if self.software_stack_changed:
                    s.append(
                        "Software environment definition has changed since last execution"
                    )

        s = "; ".join(s)
        if self.finished:
            return f"Finished (was: {s})"
        return s

    def __bool__(self):
        v = bool(
            self.updated_input
            or self.missing_output
            or self.forced
            or self.updated_input_run
            or self.noio
            or self.nooutput
            or self.pipe
            or self.service
            or self.code_changed
            or self.params_changed
            or self.software_stack_changed
            or self.input_changed
        )
        return v and not self.finished
