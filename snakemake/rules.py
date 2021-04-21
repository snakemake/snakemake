__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import re
from snakemake.path_modifier import PATH_MODIFIER_FLAG
import sys
import inspect
import sre_constants
import collections
from urllib.parse import urljoin
from pathlib import Path
from itertools import chain
from functools import partial

from snakemake.io import (
    IOFile,
    _IOFile,
    protected,
    temp,
    dynamic,
    Namedlist,
    AnnotatedString,
    contains_wildcard_constraints,
    update_wildcard_constraints,
    flag,
    get_flag_value,
)
from snakemake.io import (
    expand,
    InputFiles,
    OutputFiles,
    Wildcards,
    Params,
    Log,
    Resources,
    strip_wildcard_constraints,
)
from snakemake.io import (
    apply_wildcards,
    is_flagged,
    not_iterable,
    is_callable,
    DYNAMIC_FILL,
    ReportObject,
)
from snakemake.exceptions import (
    RuleException,
    IOFileException,
    WildcardError,
    InputFunctionException,
    WorkflowError,
    IncompleteCheckpointException,
)
from snakemake.logging import logger
from snakemake.common import Mode, lazy_property, TBDInt


class Rule:
    def __init__(self, *args, lineno=None, snakefile=None, restart_times=0):
        """
        Create a rule

        Arguments
        name -- the name of the rule
        """
        if len(args) == 2:
            name, workflow = args
            self.name = name
            self.workflow = workflow
            self.docstring = None
            self.message = None
            self._input = InputFiles()
            self._output = OutputFiles()
            self._params = Params()
            self._wildcard_constraints = dict()
            self.dependencies = dict()
            self.dynamic_output = set()
            self.dynamic_input = set()
            self.temp_output = set()
            self.protected_output = set()
            self.touch_output = set()
            self.subworkflow_input = dict()
            self.shadow_depth = None
            self.resources = None
            self.priority = 0
            self._version = None
            self._log = Log()
            self._benchmark = None
            self._conda_env = None
            self._container_img = None
            self.is_containerized = False
            self.env_modules = None
            self.group = None
            self._wildcard_names = None
            self.lineno = lineno
            self.snakefile = snakefile
            self.run_func = None
            self.shellcmd = None
            self.script = None
            self.notebook = None
            self.wrapper = None
            self.cwl = None
            self.norun = False
            self.is_handover = False
            self.is_branched = False
            self.is_checkpoint = False
            self.restart_times = 0
            self.basedir = None
            self.path_modifer = None
            self.ruleinfo = None
        elif len(args) == 1:
            other = args[0]
            self.name = other.name
            self.workflow = other.workflow
            self.docstring = other.docstring
            self.message = other.message
            self._input = InputFiles(other._input)
            self._output = OutputFiles(other._output)
            self._params = Params(other._params)
            self._wildcard_constraints = dict(other._wildcard_constraints)
            self.dependencies = dict(other.dependencies)
            self.dynamic_output = set(other.dynamic_output)
            self.dynamic_input = set(other.dynamic_input)
            self.temp_output = set(other.temp_output)
            self.protected_output = set(other.protected_output)
            self.touch_output = set(other.touch_output)
            self.subworkflow_input = dict(other.subworkflow_input)
            self.shadow_depth = other.shadow_depth
            self.resources = other.resources
            self.priority = other.priority
            self.version = other.version
            self._log = other._log
            self._benchmark = other._benchmark
            self._conda_env = other._conda_env
            self._container_img = other._container_img
            self.is_containerized = other.is_containerized
            self.env_modules = other.env_modules
            self.group = other.group
            self._wildcard_names = (
                set(other._wildcard_names)
                if other._wildcard_names is not None
                else None
            )
            self.lineno = other.lineno
            self.snakefile = other.snakefile
            self.run_func = other.run_func
            self.shellcmd = other.shellcmd
            self.script = other.script
            self.notebook = other.notebook
            self.wrapper = other.wrapper
            self.cwl = other.cwl
            self.norun = other.norun
            self.is_handover = other.is_handover
            self.is_branched = True
            self.is_checkpoint = other.is_checkpoint
            self.restart_times = other.restart_times
            self.basedir = other.basedir
            self.path_modifier = other.path_modifier
            self.ruleinfo = other.ruleinfo

    def dynamic_branch(self, wildcards, input=True):
        def get_io(rule):
            return (
                (rule.input, rule.dynamic_input)
                if input
                else (rule.output, rule.dynamic_output)
            )

        def partially_expand(f, wildcards):
            """Expand the wildcards in f from the ones present in wildcards

            This is done by replacing all wildcard delimiters by `{{` or `}}`
            that are not in `wildcards.keys()`.
            """
            # perform the partial expansion from f's string representation
            s = str(f).replace("{", "{{").replace("}", "}}")
            for key in wildcards.keys():
                s = s.replace("{{{{{}}}}}".format(key), "{{{}}}".format(key))
            # build result
            anno_s = AnnotatedString(s)
            anno_s.flags = f.flags
            return IOFile(anno_s, f.rule)

        io, dynamic_io = get_io(self)

        branch = Rule(self)
        io_, dynamic_io_ = get_io(branch)

        expansion = collections.defaultdict(list)
        for i, f in enumerate(io):
            if f in dynamic_io:
                f = partially_expand(f, wildcards)
                try:
                    for e in reversed(expand(str(f), zip, **wildcards)):
                        # need to clone the flags so intermediate
                        # dynamic remote file paths are expanded and
                        # removed appropriately
                        ioFile = IOFile(e, rule=branch)
                        ioFile.clone_flags(f)
                        expansion[i].append(ioFile)
                except KeyError:
                    return None

        # replace the dynamic files with the expanded files
        replacements = [(i, io[i], e) for i, e in reversed(list(expansion.items()))]
        for i, old, exp in replacements:
            dynamic_io_.remove(old)
            io_._insert_items(i, exp)

        if not input:
            for i, old, exp in replacements:
                if old in branch.temp_output:
                    branch.temp_output.discard(old)
                    branch.temp_output.update(exp)
                if old in branch.protected_output:
                    branch.protected_output.discard(old)
                    branch.protected_output.update(exp)
                if old in branch.touch_output:
                    branch.touch_output.discard(old)
                    branch.touch_output.update(exp)

            branch.wildcard_names.clear()
            non_dynamic_wildcards = dict(
                (name, values[0])
                for name, values in wildcards.items()
                if len(set(values)) == 1
            )
            # TODO have a look into how to concretize dependencies here
            branch._input, _, branch.dependencies = branch.expand_input(
                non_dynamic_wildcards
            )
            branch._output, _ = branch.expand_output(non_dynamic_wildcards)

            resources = branch.expand_resources(non_dynamic_wildcards, branch._input, 1)
            branch._params = branch.expand_params(
                non_dynamic_wildcards,
                branch._input,
                branch._output,
                resources,
                omit_callable=True,
            )
            branch.resources = dict(resources.items())

            branch._log = branch.expand_log(non_dynamic_wildcards)
            branch._benchmark = branch.expand_benchmark(non_dynamic_wildcards)
            branch._conda_env = branch.expand_conda_env(non_dynamic_wildcards)
            return branch, non_dynamic_wildcards
        return branch

    @property
    def is_shell(self):
        return self.shellcmd is not None

    @property
    def is_script(self):
        return self.script is not None

    @property
    def is_notebook(self):
        return self.notebook is not None

    @property
    def is_wrapper(self):
        return self.wrapper is not None

    @property
    def is_cwl(self):
        return self.cwl is not None

    @property
    def is_run(self):
        return not (
            self.is_shell
            or self.norun
            or self.is_script
            or self.is_notebook
            or self.is_wrapper
            or self.is_cwl
        )

    def check_caching(self):
        if self.name in self.workflow.cache_rules:
            if len(self.output) == 0:
                raise RuleException(
                    "Rules without output files cannot be cached.", rule=self
                )
            if len(self.output) > 1:
                prefixes = set(out.multiext_prefix for out in self.output)
                if None in prefixes or len(prefixes) > 1:
                    raise RuleException(
                        "Rules with multiple output files must define them as a single multiext() "
                        '(e.g. multiext("path/to/index", ".bwt", ".ann")). '
                        "The rationale is that multiple output files can only be unambiously resolved "
                        "if they can be distinguished by a fixed set of extensions (i.e. mime types).",
                        rule=self,
                    )
            if self.dynamic_output:
                raise RuleException(
                    "Rules with dynamic output files may not be cached.", rule=self
                )

    def has_wildcards(self):
        """
        Return True if rule contains wildcards.
        """
        return bool(self.wildcard_names)

    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, version):
        if isinstance(version, str) and "\n" in version:
            raise WorkflowError(
                "Version string may not contain line breaks.", rule=self
            )
        self._version = version

    @property
    def benchmark(self):
        return self._benchmark

    @benchmark.setter
    def benchmark(self, benchmark):
        if isinstance(benchmark, Path):
            benchmark = str(benchmark)
        if not callable(benchmark):
            benchmark = self.apply_path_modifier(benchmark, property="benchmark")
            benchmark = self._update_item_wildcard_constraints(benchmark)

        self._benchmark = IOFile(benchmark, rule=self)
        self.register_wildcards(self._benchmark.get_wildcard_names())

    @property
    def conda_env(self):
        return self._conda_env

    @conda_env.setter
    def conda_env(self, conda_env):
        self._conda_env = IOFile(conda_env, rule=self)

    @property
    def container_img(self):
        return self._container_img

    @container_img.setter
    def container_img(self, container_img):
        self._container_img = container_img

    @property
    def input(self):
        return self._input

    def set_input(self, *input, **kwinput):
        """
        Add a list of input files. Recursive lists are flattened.

        Arguments
        input -- the list of input files
        """
        for item in input:
            self._set_inoutput_item(item)
        for name, item in kwinput.items():
            self._set_inoutput_item(item, name=name)

    @property
    def output(self):
        return self._output

    @property
    def products(self):
        if self.benchmark:
            return chain(self.output, self.log, [self.benchmark])
        else:
            return chain(self.output, self.log)

    def register_wildcards(self, wildcard_names):
        if self._wildcard_names is None:
            self._wildcard_names = wildcard_names
        else:
            if self.wildcard_names != wildcard_names:
                raise SyntaxError(
                    "Not all output, log and benchmark files of "
                    "rule {} contain the same wildcards. "
                    "This is crucial though, in order to "
                    "avoid that two or more jobs write to the "
                    "same file.".format(self.name)
                )

    @property
    def wildcard_names(self):
        if self._wildcard_names is None:
            return set()
        return self._wildcard_names

    def set_output(self, *output, **kwoutput):
        """
        Add a list of output files. Recursive lists are flattened.

        After creating the output files, they are checked for duplicates.

        Arguments
        output -- the list of output files
        """
        for item in output:
            self._set_inoutput_item(item, output=True)
        for name, item in kwoutput.items():
            self._set_inoutput_item(item, output=True, name=name)

        for item in self.output:
            if self.dynamic_output and item not in self.dynamic_output:
                raise SyntaxError(
                    "A rule with dynamic output may not define any "
                    "non-dynamic output files."
                )
            self.register_wildcards(item.get_wildcard_names())
        # Check output file name list for duplicates
        self.check_output_duplicates()
        self.check_caching()

    def check_output_duplicates(self):
        """Check ``Namedlist`` for duplicate entries and raise a ``WorkflowError``
        on problems.
        """
        seen = dict()
        idx = None
        for name, value in self.output._allitems():
            if name is None:
                if idx is None:
                    idx = 0
                else:
                    idx += 1
            if value in seen:
                raise WorkflowError(
                    "Duplicate output file pattern in rule {}. First two "
                    "duplicate for entries {} and {}".format(
                        self.name, seen[value], name or idx
                    )
                )
            seen[value] = name or idx

    def apply_path_modifier(self, item, property=None):
        assert self.path_modifier is not None
        apply = partial(self.path_modifier.modify, property=property)

        assert not callable(item)
        if isinstance(item, dict):
            return {k: apply(v) for k, v in item.items()}
        elif isinstance(item, collections.abc.Iterable) and not isinstance(item, str):
            return [apply(e) for e in item]
        else:
            return apply(item)

    def update_wildcard_constraints(self):
        for i in range(len(self.output)):
            item = self.output[i]
            newitem = IOFile(
                self._update_item_wildcard_constraints(self.output[i]), rule=self
            )
            # the updated item has to have the same flags
            newitem.clone_flags(item)
            self.output[i] = newitem

    def _update_item_wildcard_constraints(self, item):
        if not (self.wildcard_constraints or self.workflow._wildcard_constraints):
            return item
        try:
            return update_wildcard_constraints(
                item, self.wildcard_constraints, self.workflow._wildcard_constraints
            )
        except ValueError as e:
            raise IOFileException(str(e), snakefile=self.snakefile, lineno=self.lineno)

    def _set_inoutput_item(self, item, output=False, name=None):
        """
        Set an item to be input or output.

        Arguments
        item     -- the item
        inoutput -- a Namedlist of either input or output items
        name     -- an optional name for the item
        """
        inoutput = self.output if output else self.input

        # Check to see if the item is a path, if so, just make it a string
        if isinstance(item, Path):
            item = str(item)
        if isinstance(item, str):
            rule_dependency = None
            if isinstance(item, _IOFile) and item.rule and item in item.rule.output:
                rule_dependency = item.rule

            item = self.apply_path_modifier(
                item, property="output" if output else "input"
            )

            # Check to see that all flags are valid
            # Note that "remote", "dynamic", and "expand" are valid for both inputs and outputs.
            if isinstance(item, AnnotatedString):
                for flag in item.flags:
                    if not output and flag in [
                        "protected",
                        "temp",
                        "temporary",
                        "directory",
                        "touch",
                        "pipe",
                    ]:
                        logger.warning(
                            "The flag '{}' used in rule {} is only valid for outputs, not inputs.".format(
                                flag, self
                            )
                        )
                    if output and flag in ["ancient"]:
                        logger.warning(
                            "The flag '{}' used in rule {} is only valid for inputs, not outputs.".format(
                                flag, self
                            )
                        )

            # add the rule to the dependencies
            if rule_dependency is not None:
                self.dependencies[item] = rule_dependency
            if output:
                item = self._update_item_wildcard_constraints(item)
            else:
                if (
                    contains_wildcard_constraints(item)
                    and self.workflow.mode != Mode.subprocess
                ):
                    logger.warning(
                        "Wildcard constraints in inputs are ignored. (rule: {})".format(
                            self
                        )
                    )
            # record rule if this is an output file output
            _item = IOFile(item, rule=self)
            if is_flagged(item, "temp"):
                if output:
                    self.temp_output.add(_item)
            if is_flagged(item, "protected"):
                if output:
                    self.protected_output.add(_item)
            if is_flagged(item, "touch"):
                if output:
                    self.touch_output.add(_item)
            if is_flagged(item, "dynamic"):
                if output:
                    self.dynamic_output.add(_item)
                else:
                    self.dynamic_input.add(_item)
            if is_flagged(item, "report"):
                report_obj = item.flags["report"]
                if report_obj.caption is not None:
                    r = ReportObject(
                        os.path.join(self.workflow.current_basedir, report_obj.caption),
                        report_obj.category,
                        report_obj.subcategory,
                        report_obj.patterns,
                        report_obj.htmlindex,
                    )
                    item.flags["report"] = r
            if is_flagged(item, "subworkflow"):
                if output:
                    raise SyntaxError("Only input files may refer to a subworkflow")
                else:
                    # record the workflow this item comes from
                    sub = item.flags["subworkflow"]
                    if _item in self.subworkflow_input:
                        other = self.subworkflow_input[_item]
                        if sub != other:
                            raise WorkflowError(
                                "The input file {} is ambiguously "
                                "associated with two subworkflows "
                                "{} and {}.".format(item, sub, other),
                                rule=self,
                            )
                    self.subworkflow_input[_item] = sub
            inoutput.append(_item)
            if name:
                inoutput._add_name(name)
        elif callable(item):
            if output:
                raise SyntaxError("Only input files can be specified as functions")
            inoutput.append(item)
            if name:
                inoutput._add_name(name)
        else:
            try:
                start = len(inoutput)
                for i in item:
                    self._set_inoutput_item(i, output=output)
                if name:
                    # if the list was named, make it accessible
                    inoutput._set_name(name, start, end=len(inoutput))
            except TypeError:
                raise SyntaxError(
                    "Input and output files have to be specified as strings or lists of strings."
                )

    @property
    def params(self):
        return self._params

    def set_params(self, *params, **kwparams):
        for item in params:
            self._set_params_item(item)
        for name, item in kwparams.items():
            self._set_params_item(item, name=name)

    def _set_params_item(self, item, name=None):
        self.params.append(item)
        if name:
            self.params._add_name(name)

    @property
    def wildcard_constraints(self):
        return self._wildcard_constraints

    def set_wildcard_constraints(self, **kwwildcard_constraints):
        self._wildcard_constraints.update(kwwildcard_constraints)

    @property
    def log(self):
        return self._log

    def set_log(self, *logs, **kwlogs):
        for item in logs:
            self._set_log_item(item)
        for name, item in kwlogs.items():
            self._set_log_item(item, name=name)

        for item in self.log:
            self.register_wildcards(item.get_wildcard_names())

    def _set_log_item(self, item, name=None):
        # Pathlib compatibility
        if isinstance(item, Path):
            item = str(item)
        if isinstance(item, str) or callable(item):
            if not callable(item):
                item = self.apply_path_modifier(item, property="log")
                item = self._update_item_wildcard_constraints(item)

            self.log.append(IOFile(item, rule=self) if isinstance(item, str) else item)
            if name:
                self.log._add_name(name)
        else:
            try:
                start = len(self.log)
                for i in item:
                    self._set_log_item(i)
                if name:
                    self.log._set_name(name, start, end=len(self.log))
            except TypeError:
                raise SyntaxError("Log files have to be specified as strings.")

    def check_wildcards(self, wildcards):
        missing_wildcards = self.wildcard_names - set(wildcards.keys())

        if missing_wildcards:
            raise RuleException(
                "Could not resolve wildcards in rule {}:\n{}".format(
                    self.name, "\n".join(self.wildcard_names)
                ),
                lineno=self.lineno,
                snakefile=self.snakefile,
            )

    def apply_input_function(
        self,
        func,
        wildcards,
        incomplete_checkpoint_func=lambda e: None,
        raw_exceptions=False,
        **aux_params
    ):
        incomplete = False
        if isinstance(func, _IOFile):
            func = func._file.callable
        elif isinstance(func, AnnotatedString):
            func = func.callable
        sig = inspect.signature(func)
        _aux_params = {k: v for k, v in aux_params.items() if k in sig.parameters}
        try:
            value = func(Wildcards(fromdict=wildcards), **_aux_params)
        except IncompleteCheckpointException as e:
            value = incomplete_checkpoint_func(e)
            incomplete = True
        except (Exception, BaseException) as e:
            if raw_exceptions:
                raise e
            else:
                raise InputFunctionException(e, rule=self, wildcards=wildcards)
        return value, incomplete

    def _apply_wildcards(
        self,
        newitems,
        olditems,
        wildcards,
        concretize=None,
        check_return_type=True,
        omit_callable=False,
        mapping=None,
        no_flattening=False,
        aux_params=None,
        apply_path_modifier=True,
        property=None,
        incomplete_checkpoint_func=lambda e: None,
        allow_unpack=True,
    ):
        if aux_params is None:
            aux_params = dict()
        for name, item in olditems._allitems():
            start = len(newitems)
            is_unpack = is_flagged(item, "unpack")
            _is_callable = is_callable(item)

            if _is_callable:
                if omit_callable:
                    continue
                item, incomplete = self.apply_input_function(
                    item,
                    wildcards,
                    incomplete_checkpoint_func=incomplete_checkpoint_func,
                    is_unpack=is_unpack,
                    **aux_params
                )
                if apply_path_modifier:
                    item = self.apply_path_modifier(item, property=property)

            if is_unpack and not incomplete:
                if not allow_unpack:
                    raise WorkflowError(
                        "unpack() is not allowed with params. "
                        "Simply return a dictionary which can be directly ."
                        "used, e.g. via {params[mykey]}."
                    )
                # Sanity checks before interpreting unpack()
                if not isinstance(item, (list, dict)):
                    raise WorkflowError(
                        "Can only use unpack() on list and dict", rule=self
                    )
                if name:
                    raise WorkflowError(
                        "Cannot combine named input file with unpack()", rule=self
                    )
                # Allow streamlined code with/without unpack
                if isinstance(item, list):
                    pairs = zip([None] * len(item), item)
                else:
                    assert isinstance(item, dict)
                    pairs = item.items()
            else:
                pairs = [(name, item)]

            for name, item in pairs:
                is_iterable = True
                if not_iterable(item) or no_flattening:
                    item = [item]
                    is_iterable = False
                for item_ in item:
                    if (
                        check_return_type
                        and not isinstance(item_, str)
                        and not isinstance(item_, Path)
                    ):
                        raise WorkflowError(
                            "Function did not return str or list " "of str.", rule=self
                        )
                    concrete = concretize(item_, wildcards, _is_callable)
                    newitems.append(concrete)
                    if mapping is not None:
                        mapping[concrete] = item_

                if name:
                    newitems._set_name(
                        name, start, end=len(newitems) if is_iterable else None
                    )
                    start = len(newitems)

    def expand_input(self, wildcards):
        def concretize_iofile(f, wildcards, is_from_callable):
            if is_from_callable:
                if isinstance(f, Path):
                    f = str(f)
                return IOFile(f, rule=self).apply_wildcards(
                    wildcards,
                    fill_missing=f in self.dynamic_input,
                    fail_dynamic=self.dynamic_output,
                )
            else:
                return f.apply_wildcards(
                    wildcards,
                    fill_missing=f in self.dynamic_input,
                    fail_dynamic=self.dynamic_output,
                )

        def handle_incomplete_checkpoint(exception):
            """If checkpoint is incomplete, target it such that it is completed
            before this rule gets executed."""
            return exception.targetfile

        input = InputFiles()
        mapping = dict()
        try:
            self._apply_wildcards(
                input,
                self.input,
                wildcards,
                concretize=concretize_iofile,
                mapping=mapping,
                incomplete_checkpoint_func=handle_incomplete_checkpoint,
                property="input",
            )
        except WildcardError as e:
            raise WildcardError(
                "Wildcards in input files cannot be " "determined from output files:",
                str(e),
                rule=self,
            )

        if self.dependencies:
            dependencies = {
                f: self.dependencies[f_]
                for f, f_ in mapping.items()
                if f_ in self.dependencies
            }
            if None in self.dependencies:
                dependencies[None] = self.dependencies[None]
        else:
            dependencies = self.dependencies

        for f in input:
            f.check()

        return input, mapping, dependencies

    def expand_params(self, wildcards, input, output, resources, omit_callable=False):
        def concretize_param(p, wildcards, is_from_callable):
            if not is_from_callable:
                if isinstance(p, str):
                    return apply_wildcards(p, wildcards)
                if isinstance(p, list):
                    return [
                        (apply_wildcards(v, wildcards) if isinstance(v, str) else v)
                        for v in p
                    ]
            return p

        params = Params()
        try:
            # When applying wildcards to params, the return type need not be
            # a string, so the check is disabled.
            self._apply_wildcards(
                params,
                self.params,
                wildcards,
                concretize=concretize_param,
                check_return_type=False,
                omit_callable=omit_callable,
                allow_unpack=False,
                no_flattening=True,
                apply_path_modifier=False,
                property="params",
                aux_params={
                    "input": input._plainstrings(),
                    "resources": resources,
                    "output": output._plainstrings(),
                    "threads": resources._cores,
                },
                incomplete_checkpoint_func=lambda e: "<incomplete checkpoint>",
            )
        except WildcardError as e:
            raise WildcardError(
                "Wildcards in params cannot be "
                "determined from output files. Note that you have "
                "to use a function to deactivate automatic wildcard expansion "
                "in params strings, e.g., `lambda wildcards: '{test}'`. Also "
                "see https://snakemake.readthedocs.io/en/stable/snakefiles/"
                "rules.html#non-file-parameters-for-rules:",
                str(e),
                rule=self,
            )
        return params

    def expand_output(self, wildcards):
        output = OutputFiles(o.apply_wildcards(wildcards) for o in self.output)
        output._take_names(self.output._get_names())
        mapping = {f: f_ for f, f_ in zip(output, self.output)}

        for f in output:
            f.check()

        # Note that we do not need to check for duplicate file names after
        # expansion as all output patterns have contain all wildcards anyway.

        return output, mapping

    def expand_log(self, wildcards):
        def concretize_logfile(f, wildcards, is_from_callable):
            if is_from_callable:
                return IOFile(f, rule=self)
            else:
                return f.apply_wildcards(
                    wildcards, fill_missing=False, fail_dynamic=self.dynamic_output
                )

        log = Log()

        try:
            self._apply_wildcards(
                log, self.log, wildcards, concretize=concretize_logfile, property="log"
            )
        except WildcardError as e:
            raise WildcardError(
                "Wildcards in log files cannot be " "determined from output files:",
                str(e),
                rule=self,
            )

        for f in log:
            f.check()

        return log

    def expand_benchmark(self, wildcards):
        try:
            benchmark = (
                self.benchmark.apply_wildcards(wildcards) if self.benchmark else None
            )
        except WildcardError as e:
            raise WildcardError(
                "Wildcards in benchmark file cannot be "
                "determined from output files:",
                str(e),
                rule=self,
            )

        if benchmark is not None:
            benchmark.check()

        return benchmark

    def expand_resources(self, wildcards, input, attempt):
        resources = dict()

        def apply(name, res, threads=None):
            if callable(res):
                aux = dict(rulename=self.name)
                if threads:
                    aux["threads"] = threads
                try:
                    try:
                        res, _ = self.apply_input_function(
                            res,
                            wildcards,
                            input=input,
                            attempt=attempt,
                            incomplete_checkpoint_func=lambda e: 0,
                            raw_exceptions=True,
                            **aux
                        )
                    except FileNotFoundError as e:
                        # Resources can depend on input files. Since expansion can happen during dryrun,
                        # where input files are not yet present, we need to skip such resources and
                        # mark them as [TBD].
                        if e.filename in input:
                            # use zero for resource if it cannot yet be determined
                            res = TBDInt(0)
                        else:
                            raise e
                except (Exception, BaseException) as e:
                    raise InputFunctionException(e, rule=self, wildcards=wildcards)

                if not isinstance(res, int) and not isinstance(res, str):
                    raise WorkflowError(
                        "Resources function did not return int or str.", rule=self
                    )
            if isinstance(res, int):
                global_res = self.workflow.global_resources.get(name, res)
                if global_res is not None:
                    res = min(global_res, res)
            return res

        threads = apply("_cores", self.resources["_cores"])
        resources["_cores"] = threads

        for name, res in self.resources.items():
            if name != "_cores":
                resources[name] = apply(name, res, threads=threads)
        resources = Resources(fromdict=resources)
        return resources

    def expand_group(self, wildcards):
        """Expand the group given wildcards."""
        if callable(self.group):
            item, _ = self.apply_input_function(self.group, wildcards)
            return item
        elif isinstance(self.group, str):
            return apply_wildcards(self.group, wildcards, dynamic_fill=DYNAMIC_FILL)
        else:
            return self.group

    def expand_conda_env(self, wildcards):
        try:
            conda_env = (
                self.conda_env.apply_wildcards(wildcards) if self.conda_env else None
            )
        except WildcardError as e:
            raise WildcardError(
                "Wildcards in conda environment file cannot be "
                "determined from output files:",
                str(e),
                rule=self,
            )

        if conda_env is not None:
            conda_env.check()

        return conda_env

    def is_producer(self, requested_output):
        """
        Returns True if this rule is a producer of the requested output.
        """
        try:
            for o in self.products:
                if o.match(requested_output):
                    return True
            return False
        except sre_constants.error as ex:
            raise IOFileException(
                "{} in wildcard statement".format(ex),
                snakefile=self.snakefile,
                lineno=self.lineno,
            )
        except ValueError as ex:
            raise IOFileException(
                "{}".format(ex), snakefile=self.snakefile, lineno=self.lineno
            )

    def get_wildcards(self, requested_output):
        """
        Return wildcard dictionary by matching regular expression
        output files to the requested concrete ones.

        Arguments
        requested_output -- a concrete filepath
        """
        if requested_output is None:
            return dict()
        bestmatchlen = 0
        bestmatch = None

        for o in self.products:
            match = o.match(requested_output)
            if match:
                l = self.get_wildcard_len(match.groupdict())
                if not bestmatch or bestmatchlen > l:
                    bestmatch = match.groupdict()
                    bestmatchlen = l
        self.check_wildcards(bestmatch)
        return bestmatch

    @staticmethod
    def get_wildcard_len(wildcards):
        """
        Return the length of the given wildcard values.

        Arguments
        wildcards -- a dict of wildcards
        """
        return sum(map(len, wildcards.values()))

    def __lt__(self, rule):
        comp = self.workflow._ruleorder.compare(self, rule)
        return comp < 0

    def __gt__(self, rule):
        comp = self.workflow._ruleorder.compare(self, rule)
        return comp > 0

    def __str__(self):
        return self.name

    def __hash__(self):
        return self.name.__hash__()

    def __eq__(self, other):
        return self.name == other.name and self.output == other.output


class Ruleorder:
    def __init__(self):
        self.order = list()

    def add(self, *rulenames):
        """
        Records the order of given rules as rule1 > rule2 > rule3, ...
        """
        self.order.append(list(rulenames))

    def compare(self, rule1, rule2):
        """
        Return whether rule2 has a higher priority than rule1.
        """
        # if rules have the same name, they have been specialized by dynamic output
        # in that case, clauses are irrelevant and have to be skipped
        if rule1.name != rule2.name:
            # try the last clause first,
            # i.e. clauses added later overwrite those before.
            for clause in reversed(self.order):
                try:
                    i = clause.index(rule1.name)
                    j = clause.index(rule2.name)
                    # rules with higher priority should have a smaller index
                    comp = j - i
                    if comp < 0:
                        comp = -1
                    elif comp > 0:
                        comp = 1
                    return comp
                except ValueError:
                    pass

        # if no ruleorder given, prefer rule without wildcards
        wildcard_cmp = rule2.has_wildcards() - rule1.has_wildcards()
        if wildcard_cmp != 0:
            return wildcard_cmp

        return 0

    def __iter__(self):
        return self.order.__iter__()


class RuleProxy:
    def __init__(self, rule):
        self.rule = rule

    @lazy_property
    def output(self):
        return self._to_iofile(self.rule.output)

    @lazy_property
    def input(self):
        return self.rule.input._stripped_constraints()

    @lazy_property
    def params(self):
        return self.rule.params._clone()

    @property
    def benchmark(self):
        return IOFile(strip_wildcard_constraints(self.rule.benchmark), rule=self.rule)

    @lazy_property
    def log(self):
        return self._to_iofile(self.rule.log)

    def _to_iofile(self, files):
        def cleanup(f):
            prefix = self.rule.workflow.default_remote_prefix
            # remove constraints and turn this into a plain string
            cleaned = strip_wildcard_constraints(f)

            modified_by = get_flag_value(f, PATH_MODIFIER_FLAG)

            if (
                self.rule.workflow.default_remote_provider is not None
                and f.startswith(prefix)
                and not is_flagged(f, "local")
            ):
                cleaned = f[len(prefix) + 1 :]
                cleaned = IOFile(cleaned, rule=self.rule)
            else:
                cleaned = IOFile(AnnotatedString(cleaned), rule=self.rule)
                cleaned.clone_remote_object(f)

            if modified_by is not None:
                cleaned.flags[PATH_MODIFIER_FLAG] = modified_by

            return cleaned

        files = Namedlist(files, custom_map=cleanup)

        return files
