__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import copy
import os
import types
import typing
from snakemake.path_modifier import PATH_MODIFIER_FLAG
import collections
from pathlib import Path
from itertools import chain
from functools import partial

try:
    import re._constants as sre_constants
except ImportError:  # python < 3.11
    import sre_constants

from snakemake_interface_executor_plugins.settings import ExecMode

from snakemake.io import (
    IOFile,
    _IOFile,
    Namedlist,
    AnnotatedString,
    contains_wildcard,
    contains_wildcard_constraints,
    get_flag_store_keys,
    is_multiext_items,
    update_wildcard_constraints,
    flag,
    get_flag_value,
    expand,
    InputFiles,
    OutputFiles,
    Wildcards,
    Params,
    Log,
    Resources,
    strip_wildcard_constraints,
    apply_wildcards,
    is_flagged,
    flag,
    is_callable,
    ReportObject,
)
from snakemake.exceptions import (
    InputOpenException,
    RuleException,
    IOFileException,
    WildcardError,
    InputFunctionException,
    WorkflowError,
    IncompleteCheckpointException,
    is_file_not_found_error,
)
from snakemake.logging import logger
from snakemake.common import (
    ON_WINDOWS,
    get_function_params,
    get_input_function_aux_params,
    mb_to_mib,
)
from snakemake.common.tbdstring import TBDString
from snakemake.resources import infer_resources
from snakemake_interface_common.utils import not_iterable, lazy_property
from snakemake_interface_common.rules import RuleInterface


_NOT_CACHED = object()


class Rule(RuleInterface):
    def __init__(self, name, workflow, lineno=None, snakefile=None):
        """
        Create a rule

        Arguments
        name -- the name of the rule
        """
        self._name = name
        self.workflow = workflow
        self.docstring = None
        self.message = None
        self._input = InputFiles()
        self._output = OutputFiles()
        self._params = Params()
        self._wildcard_constraints = dict()
        self.dependencies = dict()
        self.temp_output = set()
        self.protected_output = set()
        self.touch_output = set()
        self.shadow_depth = None
        self.resources = None
        self.priority = 0
        self._log = Log()
        self._benchmark = None
        self._conda_env = None
        self._expanded_conda_env = _NOT_CACHED
        self._container_img = None
        self.is_containerized = False
        self.env_modules = None
        self._group = None
        self._wildcard_names = None
        self._lineno = lineno
        self._snakefile = snakefile
        self.run_func = None
        self.shellcmd = None
        self.script = None
        self.notebook = None
        self.wrapper = None
        self.template_engine = None
        self.cwl = None
        self.norun = False
        self.is_handover = False
        self.is_checkpoint = False
        self._restart_times = 0
        self.basedir = None
        self.input_modifier = None
        self.output_modifier = None
        self.log_modifier = None
        self.benchmark_modifier = None
        self.ruleinfo = None
        self.module_globals = None

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def lineno(self):
        return self._lineno

    @property
    def snakefile(self):
        return self._snakefile

    @property
    def restart_times(self):
        if self.workflow.remote_execution_settings.preemptible_rules.is_preemptible(
            self.name
        ):
            return self.workflow.remote_execution_settings.preemptible_retries
        if self._restart_times is None:
            return self.workflow.execution_settings.retries
        return self._restart_times

    @restart_times.setter
    def restart_times(self, restart_times):
        self._restart_times = restart_times

    @property
    def group(self):
        if not self.workflow.non_local_exec_or_dryrun:
            return None
        else:
            overwrite_group = self.workflow.group_settings.overwrite_groups.get(
                self.name
            )
            if overwrite_group is not None:
                return overwrite_group
            return self._group

    @group.setter
    def group(self, group):
        self._group = group

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
    def is_template_engine(self):
        return self.template_engine is not None

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
        if self.workflow.cache_rules.get(self.name):
            if len(self.output) == 0:
                raise RuleException(
                    "Rules without output files cannot be cached.", rule=self
                )
            if len(self.output) > 1:
                prefixes = set(out.multiext_prefix for out in self.output)
                if None in prefixes or len(prefixes) > 1:
                    raise RuleException(
                        "Rules marked as eligible for caching that have with multiple "
                        "output files must define them as a single multiext() "
                        '(e.g. multiext("path/to/index", ".bwt", ".ann")). '
                        "The rationale is that multiple output files can only be unambiously resolved "
                        "if they can be distinguished by a fixed set of extensions (i.e. mime types).",
                        rule=self,
                    )

    def has_wildcards(self):
        """
        Return True if rule contains wildcards.
        """
        return bool(self.wildcard_names)

    @property
    def benchmark(self):
        return self._benchmark

    @benchmark.setter
    def benchmark(self, benchmark):
        if isinstance(benchmark, Path):
            benchmark = str(benchmark)
        if not callable(benchmark):
            benchmark = self.apply_path_modifier(
                benchmark, self.benchmark_modifier, property="benchmark"
            )
            benchmark = self._update_item_wildcard_constraints(benchmark)

        self._benchmark = IOFile(benchmark, rule=self)
        self.register_wildcards(self._benchmark.get_wildcard_names())

    @property
    def conda_env(self):
        return self._conda_env

    @conda_env.setter
    def conda_env(self, conda_env):
        self._conda_env = conda_env

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

        consider_ancient = self.workflow.workflow_settings.consider_ancient.get(
            self.name, frozenset()
        )
        for i, item in enumerate(input):
            if is_multiext_items(item):
                for ifile in item:
                    self._set_inoutput_item(
                        ifile,
                        name=get_flag_value(ifile, "multiext").name,
                        mark_ancient=i in consider_ancient,
                    )
            else:
                self._set_inoutput_item(
                    item,
                    mark_ancient=i in consider_ancient,
                )

        for name, item in kwinput.items():
            self._set_inoutput_item(
                item, name=name, mark_ancient=name in consider_ancient
            )

    @property
    def output(self):
        return self._output

    def products(self, include_logfiles=True):
        products = [self.output]
        if include_logfiles:
            products.append(self.log)
        if self.benchmark:
            products.append([self.benchmark])
        return chain(*products)

    def get_some_product(self):
        for product in self.products():
            return product
        return None

    def has_products(self):
        return self.get_some_product() is not None

    def register_wildcards(self, wildcard_names):
        if self._wildcard_names is None:
            self._wildcard_names = wildcard_names
        else:
            if self.wildcard_names != wildcard_names:
                raise RuleException(
                    "Not all output, log and benchmark files of "
                    "rule {} contain the same wildcards. "
                    "This is crucial though, in order to "
                    "avoid that two or more jobs write to the "
                    "same file.".format(self.name),
                    rule=self,
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
            # Named multiext have their name set under the flag (MultiextValue), if the first one is named, all of them are named.
            # Any of the output files in item can be multiext, so we do need to check all of them.
            if is_multiext_items(item):
                for ofile in item:
                    self._set_inoutput_item(
                        ofile,
                        name=get_flag_value(ofile, "multiext").name,
                        output=True,
                    )
            else:
                self._set_inoutput_item(item, output=True)
        for name, item in kwoutput.items():
            self._set_inoutput_item(item, output=True, name=name)

        for item in self.output:
            self.register_wildcards(item.get_wildcard_names())
        # Check output file name list for duplicates
        self.check_output_duplicates()
        self.check_caching()

    def check_output_duplicates(self):
        """Check ``Namedlist`` for duplicate entries and raise a ``WorkflowError``
        on problems. Does not raise if the entry is empty.
        """
        seen = dict()
        idx = None
        for name, value in self.output._allitems():
            if name is None:
                if idx is None:
                    idx = 0
                else:
                    idx += 1
            if value and value in seen:
                raise WorkflowError(
                    "Duplicate output file pattern in rule {}. First two "
                    "duplicate for entries {} and {}.".format(
                        self.name, seen[value], name or idx
                    )
                )
            seen[value] = name or idx

    def apply_path_modifier(self, item, path_modifier, property=None):
        assert path_modifier is not None
        apply = partial(path_modifier.modify, property=property)

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

            newitem = None
            if item.is_storage:
                storage_object = copy.copy(item.storage_object)
                storage_object.query = self._update_item_wildcard_constraints(
                    storage_object.query
                )
                newitem = IOFile(storage_object.local_path(), rule=self)
                newitem.clone_flags(item, skip_storage_object=True)
                newitem.flags["storage_object"] = storage_object
            else:
                newitem = IOFile(
                    self._update_item_wildcard_constraints(self.output[i]), rule=self
                )
                newitem.clone_flags(item)
            self.output[i] = newitem

    def _update_item_wildcard_constraints(self, item):
        if not (self.wildcard_constraints or self.workflow.wildcard_constraints):
            return item
        try:
            return update_wildcard_constraints(
                item, self.wildcard_constraints, self.workflow.wildcard_constraints
            )
        except ValueError as e:
            raise WorkflowError(e, snakefile=self.snakefile, lineno=self.lineno)

    def _set_inoutput_item(self, item, output=False, name=None, mark_ancient=False):
        """
        Set an item to be input or output.

        Arguments
        item     -- the item
        inoutput -- a Namedlist of either input or output items
        name     -- an optional name for the item
        """

        inoutput = self.output if output else self.input

        default_flags = (
            self.workflow.modifier.default_output_flags
            if output
            else self.workflow.modifier.default_input_flags
        )

        # Check to see if the item is a path, if so, just make it a string
        if isinstance(item, Path):
            item = str(item.as_posix())

        if isinstance(item, str):
            if ON_WINDOWS:
                if isinstance(item, (_IOFile, AnnotatedString)):
                    item = item.new_from(item.replace(os.sep, os.altsep))
                else:
                    item = item.replace(os.sep, os.altsep)

            rule_dependency = None
            if isinstance(item, _IOFile) and item.rule and item in item.rule.output:
                rule_dependency = item.rule

            if output:
                path_modifier = self.output_modifier
                property = "output"
            else:
                path_modifier = self.input_modifier
                property = "input"

            item = self.apply_path_modifier(item, path_modifier, property=property)

            item = default_flags.apply(item)

            # Check to see that all flags are valid
            # Note that "storage", and "expand" are valid for both inputs and outputs.
            if isinstance(item, AnnotatedString):
                for item_flag in item.flags:
                    if not output and item_flag in [
                        "protected",
                        "temp",
                        "nodelocal",
                        "temporary",
                        "directory",
                        "touch",
                        "pipe",
                        "service",
                        "ensure",
                        "update",
                    ]:
                        logger.warning(
                            "The flag '{}' used in rule {} is only valid for outputs, not inputs.".format(
                                item_flag, self
                            )
                        )
                    if output and item_flag in ["ancient", "before_update"]:
                        logger.warning(
                            "The flag '{}' used in rule {} is only valid for inputs, not outputs.".format(
                                item_flag, self
                            )
                        )

            if rule_dependency is not None:
                # add the rule to the dependencies
                self.dependencies[item] = rule_dependency

            if output:
                item = self._update_item_wildcard_constraints(item)
                if self.workflow.storage_settings.all_temp:
                    # mark as temp if all output files shall be marked as temp
                    item = flag(item, "temp")
            else:
                # input
                if (
                    contains_wildcard_constraints(item)
                    and self.workflow.exec_mode != ExecMode.SUBPROCESS
                ):
                    logger.warning(
                        "Wildcard constraints in inputs are ignored. (rule: {})".format(
                            self
                        )
                    )
                if mark_ancient:
                    item = flag(item, "ancient")

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
            if is_flagged(item, "report"):
                report_obj = item.flags["report"]
                if report_obj.caption is not None:
                    r = ReportObject(
                        self.workflow.current_basedir.join(report_obj.caption),
                        report_obj.category,
                        report_obj.subcategory,
                        report_obj.labels,
                        report_obj.patterns,
                        report_obj.htmlindex,
                    )
                    item.flags["report"] = r
            inoutput.append(_item)
            if name:
                inoutput._add_name(name)
        elif callable(item):
            if output:
                raise RuleException(
                    "Only input files can be specified as functions", rule=self
                )

            item = default_flags.apply(item)

            inoutput.append(item)
            if name:
                inoutput._add_name(name)
        else:
            try:
                start = len(inoutput)
                for subitem in item:
                    self._set_inoutput_item(subitem, output=output)
                if name:
                    # if the list was named, make it accessible
                    inoutput._set_name(name, start, end=len(inoutput))
            except TypeError:
                raise RuleException(
                    "Input and output files have to be specified as strings or lists of strings.",
                    rule=self,
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
                item = self.apply_path_modifier(item, self.log_modifier, property="log")
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
                raise RuleException(
                    "Log files have to be specified as strings.", rule=self
                )

    def check_wildcards(self, wildcards):
        missing_wildcards = self.wildcard_names - set(wildcards.keys())

        if missing_wildcards:
            raise RuleException(
                "Could not resolve wildcards:\n{}".format(
                    "\n".join(self.wildcard_names)
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
        groupid=None,
        **aux_params,
    ):
        if isinstance(func, _IOFile):
            func = func._file.callable
        elif isinstance(func, AnnotatedString):
            func = func.callable

        if "groupid" in get_function_params(func):
            if groupid is not None:
                aux_params["groupid"] = groupid
            else:
                # Return empty list of files and incomplete marker
                # the job will be reevaluated once groupids have been determined
                return [], True

        _aux_params = get_input_function_aux_params(func, aux_params)

        # call any callables in _aux_params
        # This way, we enable to delay the evaluation of expensive
        # aux params until they are actually needed.
        for name, value in list(_aux_params.items()):
            if callable(value):
                _aux_params[name] = value()

        wildcards_arg = Wildcards(fromdict=wildcards)

        def apply_func(func):
            incomplete = False
            try:
                value = func(wildcards_arg, **_aux_params)
                if isinstance(value, types.GeneratorType):
                    # generators should be immediately collected here,
                    # otherwise we would miss any exceptions and
                    # would have to capture them again later.
                    value = list(value)
            except IncompleteCheckpointException as e:
                value = incomplete_checkpoint_func(e)
                incomplete = True
            except InputOpenException as e:
                e.rule = self
                raise e
            except Exception as e:
                if "input" in aux_params and is_file_not_found_error(
                    e, aux_params["input"]
                ):
                    # Function evaluation can depend on input files. Since expansion can happen during dryrun,
                    # where input files are not yet present, we need to skip such cases and
                    # mark them as <TBD>.
                    value = TBDString()
                elif raw_exceptions:
                    raise e
                else:
                    raise InputFunctionException(e, rule=self, wildcards=wildcards)
            return value, incomplete

        res = func
        tries = 0
        while (callable(res) or tries == 0) and tries < 10:
            res, incomplete = apply_func(res)
            tries += 1
        if tries == 10:
            raise WorkflowError(
                "Evaluated 10 nested input functions (i.e. input functions that "
                "themselves return an input function.). More than 10 such nested "
                "evaluations are not allowed. Does the workflow accidentally return a "
                "function instead of calling it in the input function?",
                rule=self,
            )

        return res, incomplete

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
        path_modifier=None,
        property=None,
        incomplete_checkpoint_func=lambda e: None,
        allow_unpack=True,
        groupid=None,
        non_derived_items: typing.List[typing.Any] = None,
    ):
        incomplete = False
        if aux_params is None:
            aux_params = dict()
        for name, item in olditems._allitems():
            olditem = item
            start = len(newitems)
            is_unpack = is_flagged(item, "unpack")
            _is_callable = is_callable(item)
            is_derived = False

            if _is_callable:
                if omit_callable:
                    continue
                if non_derived_items is not None:
                    is_derived = self._is_deriving_function(item)
                item, incomplete = self.apply_input_function(
                    item,
                    wildcards,
                    incomplete_checkpoint_func=incomplete_checkpoint_func,
                    is_unpack=is_unpack,
                    groupid=groupid,
                    **aux_params,
                )

            if is_unpack and not incomplete:
                if not allow_unpack:
                    raise WorkflowError(
                        "unpack() is not allowed with params. "
                        "Simply return a dictionary which can be directly ."
                        "used, e.g. via {params[mykey]}.",
                        rule=self,
                    )
                # Sanity checks before interpreting unpack()
                if not isinstance(item, (list, dict)):
                    raise WorkflowError(
                        f"Can only use unpack() on list and dict, but {item} was returned.",
                        rule=self,
                    )
                if name:
                    raise WorkflowError(
                        f"Cannot combine named input file (name {name}) with unpack()",
                        rule=self,
                    )
                # Allow streamlined code with/without unpack
                if isinstance(item, list):
                    apply_results = zip(
                        [None] * len(item),
                        item,
                        [_is_callable] * len(item),
                        [is_derived] * len(item),
                    )
                else:
                    assert isinstance(item, dict)
                    apply_results = [
                        (name, item, _is_callable, is_derived)
                        for name, item in item.items()
                    ]
            else:
                apply_results = [
                    (name, item, olditem if _is_callable else None, is_derived)
                ]

            for name, item, from_callable, is_derived in apply_results:
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
                        raise InputFunctionException(
                            f"Function did not return str or iterable of str. Encountered: {item} ({type(item)})",
                            rule=self,
                            wildcards=wildcards,
                        )

                    if (
                        from_callable is not None
                        and not incomplete
                        and path_modifier is not None
                    ):
                        item_ = self.apply_path_modifier(
                            item_, path_modifier, property=property
                        )

                    concrete = concretize(item_, wildcards, from_callable)
                    newitems.append(concrete)
                    if not is_derived and non_derived_items is not None:
                        non_derived_items.append(concrete)
                    if mapping is not None:
                        mapping[concrete] = olditem

                if name:
                    newitems._set_name(
                        name, start, end=len(newitems) if is_iterable else None
                    )
                    start = len(newitems)
        return incomplete

    def expand_input(self, wildcards, groupid=None):
        def concretize_iofile(f, wildcards, from_callable):
            if from_callable is not None:
                if isinstance(f, Path):
                    f = str(f)
                iofile = IOFile(f, rule=self).apply_wildcards(wildcards)

                # inherit flags from callable
                if hasattr(from_callable, "flags"):
                    for key, value in from_callable.flags.items():
                        if key in iofile.flags:
                            continue
                        iofile.flags[key] = value

                return self.workflow.modifier.default_input_flags.apply(iofile)
            else:
                return f.apply_wildcards(wildcards)

        def handle_incomplete_checkpoint(exception):
            """If checkpoint is incomplete, target it such that it is completed
            before this rule gets executed."""
            return exception.targetfile

        input = InputFiles()
        mapping = dict()
        try:
            incomplete = self._apply_wildcards(
                input,
                self.input,
                wildcards,
                concretize=concretize_iofile,
                mapping=mapping,
                incomplete_checkpoint_func=handle_incomplete_checkpoint,
                path_modifier=self.input_modifier,
                property="input",
                groupid=groupid,
            )
        except WildcardError as e:
            raise WildcardError(
                "Wildcards in input files cannot be determined from output files:",
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

        return input, mapping, dependencies, incomplete

    @classmethod
    def _is_deriving_function(cls, func):
        if is_callable(func):
            func_params = get_function_params(func)
            return (
                "input" in func_params
                or "output" in func_params
                or "threads" in func_params
                or "resources" in func_params
            )
        return False

    def expand_params(self, wildcards, input, output, job, omit_callable=False):
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

        def handle_incomplete_checkpoint(exception):
            """If checkpoint is incomplete, target it such that it is completed
            before this rule gets executed."""
            if exception.targetfile in input:
                return TBDString()
            else:
                raise WorkflowError(
                    "Rule parameter depends on checkpoint but checkpoint output is not "
                    "defined as input file for the rule. Please add the output of the "
                    "respective checkpoint to the rule inputs. "
                    f"Input: {','.join(input)} "
                    f"Checkpoint file: {exception.targetfile}",
                    rule=self,
                )

        # We make sure that resources are only evaluated if a param function
        # actually needs them by turning them into callables and delegating their
        # evaluation to a later stage that only happens if the param function
        # requests access to resources or threads.
        resources = lambda: job.resources
        threads = lambda: job.resources._cores

        params = Params()
        non_derived_params = []
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
                property="params",
                aux_params={
                    "input": input._plainstrings(),
                    "resources": resources,
                    "output": output._plainstrings(),
                    "threads": threads,
                },
                incomplete_checkpoint_func=handle_incomplete_checkpoint,
                non_derived_items=non_derived_params,
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
        return params, non_derived_params

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
                return f.apply_wildcards(wildcards)

        log = Log()

        try:
            self._apply_wildcards(
                log,
                self.log,
                wildcards,
                concretize=concretize_logfile,
                path_modifier=self.log_modifier,
                property="log",
            )
        except WildcardError as e:
            raise WildcardError(
                "Wildcards in log files cannot be determined from output files:",
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

    def expand_resources(
        self, wildcards, input, attempt, skip_evaluation: typing.Optional[set] = None
    ):
        resources = dict()

        def apply(name, res, threads=None):
            if skip_evaluation is not None and name in skip_evaluation:
                res = TBDString()
            else:
                if isinstance(res, AnnotatedString) and res.callable:
                    res = res.callable
                if callable(res):
                    aux = dict(rulename=self.name)
                    if threads is not None:
                        aux["threads"] = threads
                    try:
                        res, _ = self.apply_input_function(
                            res,
                            wildcards,
                            input=input,
                            attempt=attempt,
                            incomplete_checkpoint_func=lambda e: 0,
                            raw_exceptions=True,
                            **aux,
                        )
                    except BaseException as e:
                        raise InputFunctionException(e, rule=self, wildcards=wildcards)

                if isinstance(res, float):
                    # round to integer
                    res = int(round(res))

                if (
                    not isinstance(res, int)
                    and not isinstance(res, str)
                    and res is not None
                ):
                    raise WorkflowError(
                        f"Resource {name} is neither int, float(would be rounded to nearest int), str, or None.",
                        rule=self,
                    )

            global_res = self.workflow.global_resources.get(name)
            if global_res is not None and res is not None:
                if not isinstance(res, TBDString) and type(res) != type(global_res):
                    global_type = (
                        "an int" if isinstance(global_res, int) else type(global_res)
                    )
                    raise WorkflowError(
                        f"Resource {name} is of type {type(res).__name__} but global resource constraint "
                        f"defines {global_type} with value {global_res}. "
                        "Resources with the same name need to have the same types (int, float, or str are allowed).",
                        rule=self,
                    )
                if isinstance(res, int):
                    res = min(global_res, res)
            return res

        threads = apply("_cores", self.resources["_cores"])
        if threads is None:
            raise WorkflowError("Threads must be given as an int", rule=self)
        if self.workflow.resource_settings.max_threads is not None and not isinstance(
            threads, TBDString
        ):
            threads = min(threads, self.workflow.resource_settings.max_threads)
        resources["_cores"] = threads

        for name, res in list(self.resources.items()):
            if name != "_cores":
                value = apply(name, res, threads=threads)

                if value is not None:
                    resources[name] = value

                    if not isinstance(value, TBDString):
                        # Infer standard resources from eventual human readable forms.
                        infer_resources(name, value, resources)
                        value = resources[name]

                    # infer additional resources
                    for mb_item, mib_item in (
                        ("mem_mb", "mem_mib"),
                        ("disk_mb", "disk_mib"),
                    ):
                        if (
                            name == mb_item
                            and mib_item not in self.resources.keys()
                            and isinstance(value, int)
                        ):
                            # infer mem_mib (memory in Mebibytes) as additional resource
                            resources[mib_item] = mb_to_mib(value)

        resources = Resources(fromdict=resources)
        return resources

    def expand_group(self, wildcards):
        """Expand the group given wildcards."""
        if callable(self.group):
            item, _ = self.apply_input_function(self.group, wildcards)
            return item
        elif isinstance(self.group, str):
            resolved = apply_wildcards(self.group, wildcards)
            if resolved != self.group:
                self.workflow.parent_groupids[resolved] = self.group
            return resolved
        else:
            return self.group

    def expand_conda_env(self, wildcards, params=None, input=None):
        if self._expanded_conda_env is not _NOT_CACHED:
            return self._expanded_conda_env

        from snakemake.common import is_local_file
        from snakemake.sourcecache import SourceFile, infer_source_file
        from snakemake.deployment.conda import (
            CondaEnvFileSpec,
            CondaEnvNameSpec,
            CondaEnvDirSpec,
            CondaEnvSpecType,
        )

        conda_env = self._conda_env
        if conda_env is not None:
            if not callable(conda_env):
                cacheable = not contains_wildcard(conda_env)
            else:
                conda_env, _ = self.apply_input_function(
                    conda_env, wildcards=wildcards, params=params, input=input
                )
                cacheable = False
                if conda_env is None:
                    return None
        else:
            self._expanded_conda_env = None
            return None

        spec_type = CondaEnvSpecType.from_spec(conda_env)

        if spec_type is CondaEnvSpecType.FILE:
            if not isinstance(conda_env, SourceFile):
                if is_local_file(conda_env) and not os.path.isabs(conda_env):
                    # Conda env file paths are considered to be relative to the directory of the Snakefile
                    # hence we adjust the path accordingly.
                    # This is not necessary in case of receiving a SourceFile.
                    conda_env = self.basedir.join(conda_env)
                else:
                    # infer source file from unmodified uri or path
                    conda_env = infer_source_file(conda_env)

            conda_env = CondaEnvFileSpec(conda_env, rule=self)
        elif spec_type is CondaEnvSpecType.NAME:
            conda_env = CondaEnvNameSpec(conda_env)
        elif spec_type is CondaEnvSpecType.DIR:
            conda_env = CondaEnvDirSpec(conda_env, rule=self)

        conda_env = conda_env.apply_wildcards(wildcards, self)
        conda_env.check()

        if cacheable:
            self._expanded_conda_env = conda_env

        return conda_env

    def expand_container_img(self, wildcards):
        """
        Expand the given container wildcards
        """
        if callable(self.container_img):
            container_url, _ = self.apply_input_function(
                self.container_img, wildcards=wildcards
            )
            return container_url

        elif isinstance(self.container_img, str):
            resolved_url = apply_wildcards(self.container_img, wildcards)
            return resolved_url

        return self.container_img

    def is_producer(self, requested_output):
        """
        Returns True if this rule is a producer of the requested output.
        """
        try:
            for o in self.products():
                if o.match(requested_output):
                    return True
            return False
        except sre_constants.error as ex:
            raise IOFileException(
                f"{ex} in wildcard statement",
                snakefile=self.snakefile,
                lineno=self.lineno,
            )
        except ValueError as ex:
            raise IOFileException(f"{ex}", snakefile=self.snakefile, lineno=self.lineno)

    def get_wildcards(self, requested_output, wildcards_dict=None):
        """
        Return wildcard dictionary by
        1. trying to format the output with the given wildcards and comparing with the requested output
        2. matching regular expression output files to the requested concrete ones.

        Arguments
        requested_output -- a concrete filepath
        """
        if requested_output is None:
            return dict()

        # first try to match the output with the given wildcards
        if wildcards_dict is not None:
            if self.wildcard_names <= wildcards_dict.keys():
                wildcards_dict = {
                    wildcard: value
                    for wildcard, value in wildcards_dict.items()
                    if wildcard in self.wildcard_names
                }
                for o in self.products():
                    try:
                        applied = o.apply_wildcards(wildcards_dict)
                        # if the output formatted with the wildcards matches the requested output,
                        if applied == requested_output:
                            # we check whether the wildcards satisfy the constraints
                            constraints = o.wildcard_constraints()

                            def check_constraint(wildcard, value):
                                constraint = constraints.get(wildcard)
                                return constraint is None or constraint.match(value)

                            if all(
                                check_constraint(name, value)
                                for name, value in wildcards_dict.items()
                            ):
                                # and then just return the given wildcards_dict limited to the wildcards that are actually used
                                return wildcards_dict
                    except WildcardError:
                        continue

        bestmatchlen = 0
        bestmatch = None

        for o in self.products():
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
        if isinstance(other, Rule):
            return self.name == other.name and self.output == other.output
        else:
            return False


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
        def modify_callable(item):
            if is_callable(item):
                # For callables ensure that the rule's original path modifier is applied as well.

                def inner(wildcards):
                    return self.rule.apply_path_modifier(
                        item(wildcards), self.rule.input_modifier, property="input"
                    )

                return inner
            else:
                # For strings, the path modifier has been already applied.
                return item

        return InputFiles(
            toclone=self.rule.input, strip_constraints=True, custom_map=modify_callable
        )

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
            prefix = self.rule.workflow.storage_settings.default_storage_prefix
            # remove constraints and turn this into a plain string
            cleaned = strip_wildcard_constraints(f)

            modified_by = get_flag_value(f, PATH_MODIFIER_FLAG)

            if (
                self.rule.workflow.storage_settings.default_storage_provider is not None
                and f.startswith(prefix)
                and not is_flagged(f, "local")
            ):
                start = len(prefix) + 1 if prefix else 0
                cleaned = f.storage_object.query[start:]
                cleaned = IOFile(cleaned, rule=self.rule)
            else:
                cleaned = IOFile(AnnotatedString(cleaned), rule=self.rule)
                cleaned.clone_storage_object(f)

            if modified_by is not None:
                cleaned.flags[PATH_MODIFIER_FLAG] = modified_by

            return cleaned

        files = Namedlist(files, custom_map=cleanup)

        return files
