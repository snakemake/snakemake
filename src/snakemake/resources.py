from __future__ import annotations
from collections import defaultdict
import copy
import functools as ft
import itertools as it
import operator as op
import os
import re
import shutil
import tempfile
import math
from typing import (
    TYPE_CHECKING,
    Any,
    Awaitable,
    Callable,
    Collection,
    Dict,
    Final,
    Iterable,
    Iterator,
    List,
    Literal,
    Mapping,
    Optional,
    Self,
    Sequence,
    Tuple,
    TypeAlias,
    TypeVar,
    cast,
)

import humanfriendly
from humanfriendly import InvalidTimespan, InvalidSize, parse_size, parse_timespan
from snakemake_interface_storage_plugins.exceptions import FileOrDirectoryNotFoundError

from snakemake.common import (
    get_input_function_aux_params,
    mb_to_mib,
    mib_to_mb,
)
from snakemake.exceptions import (
    NestedCoroutineError,
    ResourceConstraintError,
    ResourceDuplicationError,
    ResourceError,
    ResourceInsufficiencyError,
    ResourceScopesException,
    ResourceConversionError,
    ResourceValidationError,
    WorkflowError,
    is_file_not_found_error,
)
from snakemake.common.tbdstring import TBDString
from snakemake.io import AnnotatedString

if TYPE_CHECKING:
    from snakemake.jobs import Job
    from snakemake.io import Wildcards
    from snakemake.settings.types import ValidResource


class GroupResources:
    @classmethod
    def basic_layered(
        cls,
        toposorted_jobs: List[List[Job]],
        constraints: Resources,
        run_local: bool,
        additive_resources: Optional[List[str]] = None,
        sortby: Optional[List[str]] = None,
    ) -> Mapping[str, str | int]:
        """Basic implementation of group job resources calculation

        Each toposort level is individually sorted into a series of layers, where each
        layer fits within the constraints. Resource constraints represent a "width" into
        which the layer must fit. For instance, with a mem_mb constraint of 5G, all the
        jobs in a single layer must together not consume more than 5G of memory. Any
        jobs that would exceed this constraint are pushed into a new layer. The overall
        width for the entire group job is equal to the width of the widest layer.

        Additive resources (by default, "runtime") represent the "height" of the layer.
        They are not directly constrained, but their value will be determined by the
        sorting of jobs based on other constraints. Each layer's height is equal to the
        height of its tallest job. For instance, a layer containing a 3hr job will have
        a runtime height of 3 hr. The total height of the entire group job will be the
        sum of the heights of all the layers.

        Note that both height and width are multidimensial, so layer widths will be
        calculated with respect to every constraint created by the user.

        In this implementation, there is no mixing of layers, which may lead to "voids".
        For instance, a layer containing a tall, skinny job of 3hr length and 1G mem
        combined with a short, fat job of 10min length and 20G memory would have a 2hr
        50min period where 19G of memory are not used. In practice, this void will be
        filled by the actual snakemake subprocess, which performs real-time scheduling
        of jobs as resources become available.  But it may lead to overestimation of
        resource requirements.

        To help mitigate against voids, this implementation sorts the jobs within a
        toposort level before assignment to layers. Jobs are first sorted by their
        overall width relative to the available constraints. So the fattest jobs will
        grouped together on the same layer. Jobs are then sorted by the resources
        specified in ``sortby``, by default "runtime". So jobs of similar length will be
        grouped on the same layer.

        Users can help mitigate against voids by grouping jobs of similar resource
        dimensions.  Eclectic groups of various runtimes and resource consumptions will
        not be estimated as efficiently as groups of homogeneous consumptions.

        Parameters
        ----------
        toposorted_jobs : list of lists of jobs
            Jobs sorted into toposort levels: the jobs in each level only depend on jobs
            in previous levels.
        constraints : dict of str -> int
            Upper limit of resource allowed. Resources without constraints will be
            treated as infinite
        run_local : bool
            True if the group is being run in the local process, rather than being
            submitted. Relevant for Pipe groups and Service groups
        additive_resources : list of str, optional
            Resources to be treated as the "height" of each layer, i.e. to be summed
            across layers.
        sortby : list of str, optional
            Resources by which to sort jobs prior to layer assignment.

        Returns
        -------
        Dict of str -> int,str
            Total resource requirements of the group job

        Raises
        ------
        WorkflowError
            Raised if an individual job requires more resources than the constraints
            allow (chiefly relevant for pipe groups)
        """
        additive_resources = (
            additive_resources if additive_resources is not None else ["runtime"]
        )
        sortby = sortby if sortby is not None else ["runtime"]

        blocks: List[Dict[str, str | int]] = []
        # iterate over siblings that can be executed in parallel
        for siblings in toposorted_jobs:
            # Total resource requirements for this toposort layer
            block_resources: Dict[str, str | int] = {}

            job_resources: List[Dict[str, str | int]] = []
            pipe_resources: Dict[str, List[Dict[str, str | int]]] = defaultdict(list)
            for job in siblings:
                # Get resources, filtering out FileNotFoundErrors. List items will
                # be job resources objects with (resource: value)
                # [
                #   { "runtime": 5, "threads": 2, "tmpdir": "/tmp" },
                #   { "runtime": 15, "tmpdir": "/tmp"},
                #   ...
                # ]
                # Pipe jobs and regular jobs are put in separate lists.
                try:
                    # Remove any TBDStrings from values. These will typically arise
                    # here because the default mem_mb and disk_mb are based off of
                    # input file size, and intermediate files in the group are not yet
                    # generated. Thus rules consuming such files must explicitly
                    # specify their resources
                    res = {
                        k: res
                        for k, res in job.resources.items()
                        if not isinstance(res, TBDString) and k not in SizedResources
                    }
                    if job.pipe_group:
                        pipe_resources[job.pipe_group].append(res)
                    else:
                        job_resources.append(res)
                except (FileNotFoundError, FileOrDirectoryNotFoundError):
                    # Skip job if resource evaluation leads to a file not found error.
                    # This will be caused by an inner job, which needs files created by
                    # the same group. All we can do is to ignore such jobs for now.
                    continue

            # Jobs in pipe groups must be run simultaneously, so we merge all the
            # resources of each pipe group into one big "job". Resources are combined
            # as "intralayer", so additives (like runtime) are maxed, and the rest are
            # summed
            for pipe in pipe_resources.values():
                job_resources.append(
                    cls._merge_resource_dict(
                        pipe,
                        methods={res: max for res in additive_resources},
                        default_method=sum,
                    )
                )

            # Set of resource types requested in at least one job
            resource_types = list(set(it.chain(*job_resources)))
            int_resources: Dict[str, List[int]] = {}
            # Sort all integer resources in job_resources into int_resources. Resources
            # defined as a string are placed immediately into block_resources.
            for res in resource_types:
                if res == "_nodes":
                    continue

                values = [resources.get(res, 0) for resources in job_resources]

                if cls._is_string_resource(res, values):
                    block_resources[res] = cast(str, values[0])
                else:
                    int_resources[res] = cast(List[int], values)

            # Collect values from global_resources to use as constraints.
            sorted_constraints = {
                name: constraints.get(name).value for name in int_resources
            }

            # For now, we are unable to handle a constraint on runtime, so ignore.
            # Jobs requesting too much runtime will still get flagged by the
            # scheduler
            for res in additive_resources:
                if res in sorted_constraints:
                    sorted_constraints[res] = None

            # Get layers
            try:
                layers = cls._get_layers(
                    int_resources, list(sorted_constraints.values()), sortby
                )
            except ResourceError as err:
                raise ResourceInsufficiencyError(
                    additive_resources, err.args[0]
                ) from err

            # Merge jobs within layers
            intralayer_merge_methods = [
                max if res in additive_resources else sum for res in int_resources
            ]
            merged = [
                cls._merge_resource_layer(layer, intralayer_merge_methods)
                for layer in layers
            ]

            # Combine layers
            interlayer_merge_methods = [
                sum if res in additive_resources else max for res in int_resources
            ]
            combined = cls._merge_resource_layer(merged, interlayer_merge_methods)

            # Reassign the combined values from each layer to their resource names
            block_resources.update(dict(zip(int_resources, combined)))

            blocks.append(block_resources)

        if run_local:
            final = cls._merge_resource_dict(blocks, default_method=sum)
        else:
            final = cls._merge_resource_dict(
                blocks,
                default_method=max,
                methods={res: sum for res in additive_resources},
            )

        final["_nodes"] = 1

        for resource in SizedResources:
            if f"{resource}_mb" in final:
                final[resource] = Resource(
                    resource, final[f"{resource}_mb"]
                ).format_human_friendly()

        return final

    @classmethod
    def _is_string_resource(cls, name: str, values: List[str | int]):
        # If any one of the values provided for a resource is not an int, we
        # can't process it in any way. So we constrain all such resource to be
        # the same
        if all([isinstance(val, int) for val in values]):
            return False
        else:
            unique = set(values)
            if len(unique) > 1:
                raise WorkflowError(
                    "Resource {name} is a string but not all jobs in group require the "
                    "same value. Observed values: {values}.".format(
                        name=name, values=unique
                    )
                )
            return True

    @classmethod
    def _merge_resource_dict(
        cls,
        resources: List[Dict[str, str | int]],
        skip: Optional[List[str]] = None,
        methods: Optional[Dict[str, Callable[[List[int]], int]]] = None,
        default_method: Callable[[List[int]], int] = max,
    ):
        skip = skip or []
        methods = methods or {}
        grouped: Dict[str, List[str | int]] = {}
        for job in resources:
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

        ret: Dict[str, int | str] = {}
        for res, values in grouped.items():
            if res in skip:
                continue

            if cls._is_string_resource(res, values):
                ret[res] = cast(str, values[0])
            elif res in methods:
                ret[res] = methods[res](cast(List[int], values))
            else:
                ret[res] = default_method(cast(List[int], values))
        return ret

    @classmethod
    def _merge_resource_layer(
        cls,
        resources: Sequence[Tuple[int, ...]],
        methods: Sequence[Callable[[Sequence[int]], int]],
    ):
        """
        Sum or max across all resource types within a layer, similar to
        summing along axis 0 in numpy:
        [
          ( 3 ^ , 4 ^ , 1 ),
          ( 2 | , 1 | , 6 ),
          ( 1 | , 4 | , 0 ),
        ]
        The method for each column is specified in methods, which should be an array
        with one index per column
        """
        return [method(r) for method, r in zip(methods, zip(*resources))]

    @staticmethod
    def _check_constraint(
        resources: List[Tuple[int, ...]], constraints: List[int | None]
    ):
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

    @classmethod
    def _get_layers(
        cls,
        resources: Dict[str, List[int]],
        constraints: List[int | None],
        sortby: Optional[List[str]] = None,
    ):
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
        def _proportion(group: Tuple[int, ...]):
            return [r / c if c else 0 for r, c in zip(group, constraints)]

        # Return the highest _proportion item in the list
        def _highest_proportion(group: Tuple[int, ...]):
            return max(_proportion(group))

        # Rows should always have at least one empty row to ensure space for insertion.
        rows: List[List[Tuple[int, ...]]] = [[]]

        # By zipping, we combine the vals into tuples based on job, 1 tuple per
        # job: [ (val1, 1_val1, 2_val1), ...]. In each tuple, the resources
        # will remain in the same order as the original dict, so their identity
        # can be extracted later.
        resource_groups: zip[Tuple[int, ...]] = zip(*resources.values())

        # Sort by _proportion highest to lowest
        pre_sorted = sorted(resource_groups, key=_highest_proportion, reverse=True)

        # If a key is provided (e.g. runtime), we sort again by that key from
        # highest to lowest
        for res in sortby or []:
            if res in resources:
                # Find the position of the key in the job tuple
                i = list(resources).index(res)
                pre_sorted = sorted(pre_sorted, key=op.itemgetter(i), reverse=True)

        for group in pre_sorted:
            appended = False

            # Check each row for space, starting with the first.
            for row in rows:
                if not appended and cls._check_constraint(row + [group], constraints):
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
                too_high: List[Tuple[str, int, int | None]] = []
                for i, val in enumerate(_proportion(group)):
                    if val > 1:
                        too_high.append((list(resources)[i], group[i], constraints[i]))

                error_text = [
                    f"\t{res}: {amount}/{constraint}"
                    for res, amount, constraint in too_high
                ]
                raise ResourceError("\n".join(error_text))

        # Remove final empty row. (The above loop ends each cycle by ensuring
        # there's an empty row)
        rows.pop()
        return rows


def mb_to_str(size: int) -> str:
    return humanfriendly.format_size(size * (10**6))


SizedResources = {"mem", "disk"}
TimeResources = {"runtime"}
HumanFriendlyResources = TimeResources | SizedResources


_T = TypeVar("_T")


class Resource:
    """Standardized representation of a resource name and value.

    In addition to ``int`` and ``str``, the standard resource types, ``Resource`` can
    store ``None``, signifying an unset resource, and a callable, signifying a resource
    to be calculated in the context of a rule. Via the ``.value`` property, ``Resource``
    ensures a resource has been computed before its value is retrieved.

    The class additionally handles input validation and standardization.
    ``float``s get rounded into ``int``. If the resource name corresponds to a size
    resource ("disk" and "mem") or a timespan resource ("runtime"), strings values
    will be interpreted as a human readable resource and converted to an integer of
    appropriate denomination (Mb for size, hr for timespan).

    Arguments
    =========
    name: str
        Name of the resource (note that "disk", "mem", and "runtime" has special
        behaviour, as described above).
    value: int, str, Callable, None
        The concrete or unevaluated value

    Raises
    ======
    ResourceValidationError:
        if the resource is not a ``str``, ``int``, ``float``, ``None``, or a callable.
    WorkflowError
        if the resource is of the human-readable group but cannot be parsed
    """

    _evaluator: Final[Callable[..., ValidResource] | None]

    def __init__(self, name: str, value: ValidResource, raw: int | str | None = None):
        if not (
            isinstance(value, (str, int, float)) or callable(value) or value is None
        ):
            msg = (
                f"Resource '{name}' assigned invalid value {value!r}. Must be str, "
                "int, float, or callable (function)."
            )
            raise ResourceValidationError(msg)
        if isinstance(value, float):
            value = round(value)
        self.name = name
        self._value = self.parse_human_friendly(name, value)

        if raw is None and isinstance(self._value, (int, str)):
            self.raw = self._value
        else:
            self.raw = raw

        if callable(self._value):
            self._evaluator = self._value
        elif isinstance(self._value, AnnotatedString) and self._value.is_callable():
            self._evaluator = cast(Callable[..., ValidResource], self._value.callable)
        else:
            self._evaluator = None

    def __repr__(self):
        if self.is_evaluable():
            value = "function(...)"
        elif isinstance(self._value, str):
            value = f'"{self._value}"'
        else:
            value = self._value
        return f'Resource("{self.name}", {value})'

    def __eq__(self, other: object):
        val = self._value
        if isinstance(other, Resource):
            return val == other._value
        else:
            return val == other

    @property
    def value(self) -> str | int | None:
        """Safely return concrete resource values.

        Errors if the resource has not yet been evaluated.
        """
        if self.is_evaluable():
            raise ValueError(
                "Resource must be evaluated within the context of a rule to "
                "compute its value"
            )
        return cast(str | int | None, self._value)

    def is_evaluable(self) -> bool:
        """Indicates if a resource has yet to be evaluated."""
        return self._evaluator is not None

    def evaluate(self, *args: Any, **kwargs: Any):
        """Evaluate the resource with given args.

        This method ensures the returned resource is correctly returned in a
        new ``Resource`` class, validated, and standardized. Resources should thus
        only be evaluated via this method.
        """
        if self._evaluator is None:
            return self
        kept_args = get_input_function_aux_params(self._evaluator, kwargs)
        try:
            return self.__class__(
                self.name, self._evaluator(*args, **kept_args), raw=self.raw
            )
        except Exception as err:
            if is_file_not_found_error(err, kwargs.get("input", tuple())):
                return self.__class__(self.name, TBDString(), raw=self.raw)
            raise err

    def constrain(self, other: Resource | int | None):
        """Use ``other`` as the maximum value for ``Self``, but only if both are integers.

        Returns self unaltered if either resource is ``None`` or if both are the same
        non-integer type.

        Raises
        ======
        ResourceConstraintError:
            if ``other`` has a different type from ``Self``, or if either are evaluable
            resources.
        """
        if isinstance(self._value, TBDString):
            return self

        other_val = other.value if isinstance(other, Resource) else other

        if self._value is None or other_val is None:
            return self

        self_val = self.value
        if not isinstance(self_val, type(other_val)):
            msg = (
                f"{self} of type '{type(self_val).__name__}' was "
                f"constrained by mismatched type '{type(other_val).__name__}' with "
                f"value '{other_val}'."
            )
            raise ResourceConstraintError(msg)
        if isinstance(self_val, int):
            assert isinstance(other_val, int)
            return Resource(self.name, min(other_val, self_val))
        return self

    def from_mb_to_mib(self):
        """Convert the resource to mebibytes.

        The name of the resource is updated to have a _mib suffix.

        The method does not track the current units of the resource, so will not prevent
        mebibytes from being "converted" into mebibytes multiple times.

        Un-evaluated resources can be converted. This is to support resource
        normalization within ``Resources`` class: units can be converted before
        the resource is evaluated in the context of a rule.

        Raises
        ======
        ResourceTypeError:
            if conversion is attempted on a str resource (or a callable that returns a
            str/None).
        """
        return self.__class__(
            f"{self.name.removesuffix('_mb')}_mib",
            self._convert_units(mb_to_mib),
        )

    def format_human_friendly(self) -> str | TBDString:
        """Convert the resource into a human friendly string.

        The units are ASSUMED to be in MB. Resource does not internally track units, and
        this method will work (incorrectly) if called an a unit suffixed with _mib.

        An _mb suffix will be stripped, if found (_mib suffixes are left untouched to
        avoid hiding the effects of an uninintended conversion.)

        Un-evaluated resources can be converted. This is to support resource
        normalization within ``Resources`` class: units can be converted before
        the resource is evaluated in the context of a rule.

        Raises
        ======
        ResourceTypeError:
            if conversion is attempted on a str resource (or a callable that returns a
            str/None).
        """
        if self.is_evaluable():
            # Use builtin error type as this should not happen
            raise TypeError(
                f"{self.name} formatted as human-friendly before evaluation"
            )
        return self._convert_units(mb_to_str)  # type: ignore (already checked not evaluable)

    def _convert_units(self, converter: Callable[[int], _T]):
        """Perform a unit conversion on an integer or callable."""
        if isinstance(self._value, TBDString):
            return self._value
        elif self.is_evaluable():
            return self._wrap_evaluator(converter)
        elif not isinstance(self._value, int):
            raise ResourceConversionError.format(self.name, self._value)
        return converter(self._value)

    def _wrap_evaluator(
        self,
        converter: Callable[[int], _T],
    ):
        assert self._evaluator is not None, "self._evaluator is None"

        @ft.wraps(self._evaluator)
        def inner(*args: Any, **kwargs: Any):
            result = self._evaluator(*args, **kwargs)  # type: ignore (self._evaluator should never change)
            if not isinstance(result, (int, float)):
                raise ResourceConversionError.format_evaluated(self.name, result)

            return converter(round(result))

        return inner

    def standardize_size(self):
        """Standardize the representation of sized resources.

        If the resource is evaluable and ends in _mib or _mb, the evaluable will be
        wrapped to ensure it returns an int or float, and mib will be converted to mb as
        appropriate.

        Byte resources are always stored with a _mb suffix.

        Raises
        ======
        ResourceTypeError:
            if conversion is attempted on a str resource (or a callable that returns a
            str/None).
        """
        if self.name.endswith("_mb"):
            # using _convert_units with the identity lets us validate that the resource
            # was set with an int/float or an evaluable returning an int/float
            return self.__class__(
                self.name, self._convert_units(lambda x: x), raw=self.raw
            )
        if self.name.endswith("_mib"):
            return self.__class__(
                self.name[:-4] + "_mb", self._convert_units(mib_to_mb), raw=self.raw
            )
        if self._evaluator is not None:
            # Wrap the evaluable so that it can be interpreted in the context of the
            # original unprefixed resource and potentially return a humanreadable value
            @ft.wraps(self._evaluator)
            def wrapper(*args: Any, **kwargs: Any):
                return self.evaluate(*args, **kwargs).value

            value = wrapper
        else:
            value = self._value

        return self.__class__(self.name + "_mb", value)

    def with_name(self, name: str):
        """Update the name of the resource without changing the value."""
        cp = copy.copy(self)
        cp.name = name
        return cp

    @classmethod
    def from_cli_expression(
        cls, name: str, value: str, *, with_threads_arg: bool = True
    ):
        """Create a new evaluable resource based on a python expression.

        Threads can optionally be included in the expression environment

        Arguments
        =========
        name:
            name of the resource
        value:
            python expression to be evaluated
        with_threads_arg: boolean
            If True, include ``threads`` as an argument in the returned function
        """

        def threads_evaluator(
            wildcards: Wildcards,
            input: Any,
            attempt: int,
            threads: int,
            rulename: str,
            async_run: Callable[[Awaitable[_T]], _T],
        ):
            return cls.cli_evaluator(
                name,
                value,
                input=input,
                attempt=attempt,
                async_run=async_run,
                threads=threads,
            )

        def nonthreads_evaluator(
            wildcards: Wildcards,
            input: Any,
            attempt: int,
            rulename: str,
            async_run: Callable[[Awaitable[_T]], _T],
        ):
            return cls.cli_evaluator(
                name, value, input=input, attempt=attempt, async_run=async_run
            )

        if with_threads_arg:
            return Resource(name, threads_evaluator, raw=value)

        return Resource(name, nonthreads_evaluator, raw=value)

    @staticmethod
    def cli_evaluator(
        name: str,
        val: str,
        *,
        input: Any,
        attempt: int,
        threads: int | None = None,
        async_run: Callable[[Awaitable[_T]], _T],
    ):
        """Evaluate a python expression to determine resource value.

        Intended to be used from ``from_cli_expression``.
        """
        # Expand env variables
        val = os.path.expanduser(os.path.expandvars(val))
        try:
            value = eval(
                val,
                {
                    "input": input,
                    "attempt": attempt,
                    "system_tmpdir": tempfile.gettempdir(),
                    "shutil": shutil,
                    "async_run": async_run,
                    **({"threads": threads} if threads is not None else {}),
                },
            )
        # Triggers for string arguments like n1-standard-4
        except (NameError, SyntaxError):
            return val
        except NestedCoroutineError:
            raise
        except Exception as e:
            # We need to repeat this logic from self.evaluate in order to support the
            # custom error message below.
            if is_file_not_found_error(e, input):
                return TBDString()

            if is_ordinary_string(val):
                return val

            raise WorkflowError(
                f"Failed to evaluate resources value '{val}'.",
                "When interpreted as a python expression, the following error was "
                "given:",
                "",
                e,
                "",
                "Additional quoting may be needed to force interpretation as a "
                "string. E.g.: --default-resources "
                "\"tmpdir='/home/user/tmp'\" or "
                "--set-resources \"somerule:someresource='--nice=100'\". "
                "This also holds for setting resources inside of a profile, where "
                "you might have to enclose them in single and double quotes, "
                "i.e. someresource: \"'--nice=100'\".",
            )
        return value

    @staticmethod
    def parse_human_friendly(name: str, value: ValidResource) -> ValidResource:
        """Parse a resource as human friendly if allowable for the resource name.

        If a non-string value is given, it is returned without modification.

        Only human-friendly resources (e.g. "mem", "disk", "runtime") are parsed, others
        are returned without modification.

        Raises
        ======
        WorkflowError:
            if a valid human-friendly resource is given but it cannot be parsed.
        """
        if isinstance(value, str) and not isinstance(value, TBDString):
            stripped = value.strip("'\"")
            err_msg = (
                f"Resource '{name}' with value {value!r} could not be parsed as "
                "{unit}"
            )
            if name in SizedResources:
                try:
                    return max(int(math.ceil(parse_size(stripped) / 1e6)), 1)
                except InvalidSize as err:
                    raise WorkflowError(err_msg.format(unit="size in MB")) from err
            elif name in TimeResources:
                try:
                    return max(int(round(parse_timespan(stripped) / 60)), 1)
                except InvalidTimespan as err:
                    raise WorkflowError(err_msg.format(unit="minutes")) from err
        return value


class Resources(Mapping[str, Resource]):
    """Standardized container of resources.

    Implements and enforces standardization of size resources (e.g. disk, disk_mb,
    mem_mib, etc). ``Resources`` cannot be initialized with both suffixed and
    unsuffixed versions of the same size resource (e.g. disk and disk_mb). If a suffixed
    resource is provided, it will be converted into the unsuffixed version with
    appropriate unit conversion (e.g. disk_mib will be saved within ``Resources`` as
    disk and its value will be unit-converted).

    Additionally implements the default resource feature of snakemake

    Initialization should be performed either with ``parse`` for CLI-style string
    assignments, or ``from_mapping``.

    Raises
    ======
    ResourceDuplicationError
        if multiple variants of the same sized resource are provided (e.g. mem, mem_mb)
    ResourceTypeError:
        if a human readable resource with an incorrect type is given (e.g. a string to a
        suffixed resource)
    """

    DEFAULTS = {
        "full": {
            "mem": "min(max(2*input.size_mb, 1000), 8000)",
            "disk": "max(2*input.size_mb, 1000) if input else 50000",
            "tmpdir": "system_tmpdir",
        },
        "bare": {"tmpdir": "system_tmpdir"},
    }

    def __init__(self, mapping: Dict[str, Resource] | None = None):
        if mapping is None:
            self._data: Dict[str, Resource] = {}
            return
        self._data = mapping
        for res in SizedResources:
            self._normalize_sizes(res)

    def _normalize_sizes(self, resource: str):
        found: set[str] = set()
        for suffix in {"", "_mb", "_mib"}:
            if resource + suffix in self._data:
                found.add(resource + suffix)
        if not found:
            return
        if len(found) > 1:
            raise ResourceDuplicationError(list(found))
        name = found.pop()
        standardized = self._data.pop(name).standardize_size()
        self._data[standardized.name] = standardized

    def __repr__(self):
        return str(self._data)

    def __iter__(self):
        return iter(self._data.keys())

    def __len__(self):
        return len(self._data)

    def __getitem__(self, ix: str, /):
        if not isinstance(ix, str):
            raise TypeError(f"{ix} must be a str")
        return self._data[ix]

    def __setitem__(self, ix: str, val: ValidResource | Resource | None, /):
        if isinstance(val, Resource):
            self._data[ix] = val
            return
        self._data[ix] = Resource(ix, val)

    @property
    def args(self):
        return [
            f"{name}={value.raw}"
            for name, value in self._data.items()
            if value.raw is not None
        ]

    @classmethod
    def decode_arg(cls, arg: str) -> Tuple[str, str]:
        try:
            return cast(Tuple[str, str], tuple(arg.split("=", maxsplit=1)))
        except ValueError:
            raise ValueError("Resources have to be defined as name=value pairs.")

    @classmethod
    def default(cls, mode: Literal["bare"] | Literal["full"] | None = None):
        """Return an object loaded with a standard default resource set.

        Arguments
        ========
        mode: "bare" or "full"
            full presets "tmpdir", "mem_mb", and "disk_mb". bare sets only "tmpdir"
        """
        try:
            default_resources = cls.DEFAULTS[mode or "full"]
        except KeyError:
            raise ValueError(f"Unexpected mode for DefaultResources: {mode}")

        result = {
            res: Resource.from_cli_expression(res, val, with_threads_arg=True)
            for res, val in default_resources.items()
        }

        return cls(result)

    @classmethod
    def parser_factory(
        cls,
        *,
        allow_expressions: bool = False,
        defaults: None | Literal["bare"] | Literal["full"] = None,
    ) -> Callable[[List[str]], Self]:
        """Return a parsing function with preset keyword arguments.

        Intended for use with argparse.
        """
        return lambda exprs: cls.parse(
            exprs,
            allow_expressions=allow_expressions,
            defaults=defaults,
        )

    @classmethod
    def parse(
        cls,
        exprs: List[str],
        *,
        allow_expressions: bool,
        defaults: None | Literal["bare"] | Literal["full"] = None,
    ):
        """Parse a series of CLI-style string resource assignments.

        Arguments
        =========
        exprs: list of str
            Each expr is in the form RESOURCE=VALUE, where RESOURCE is a valid python
            identifier, and VALUE is a number, a quoted string, or a valid python
            expression.
        allow_expressions:
            Allows use of quoted strings and python expressions in the list of exprs.
            Otherwise, an error is raised if these are found.
        defaults: "bare" or "full"
            Initializes the returned resource object with the given default set, if
            provided. Provided exprs will override default resources.
        """
        valid = re.compile(r"[a-zA-Z_]\w*$")
        unparsed = dict([cls.decode_arg(expr) for expr in exprs])

        result: Dict[str, Resource] = {}
        for res, val in unparsed.items():
            if res == "_cores":
                raise WorkflowError(
                    "Resource _cores is already defined internally. Use a different "
                    "name."
                )
            if not valid.match(res):
                raise WorkflowError(
                    "Resource definition must start with a valid identifier, but found "
                    "{}.".format(res)
                )

            # check if resource is parsable as human friendly (given the correct
            # name and formatted value). If it is, we return the parsed value to
            # save a step later.
            try:
                val = Resource.parse_human_friendly(res, val)
            except WorkflowError:
                pass

            try:
                val = int(val)
            except ValueError:
                pass
            else:
                if val < 0:
                    msg = (
                        f"Resource '{res}' was given value '{val}', which is not a "
                        f"positive integer"
                    )
                    raise WorkflowError(msg)
                result[res] = Resource(res, val)
                continue

            if allow_expressions:
                result[res] = Resource.from_cli_expression(
                    res, val, with_threads_arg=True
                )
                continue

            raise WorkflowError(
                f"Resource '{res}' was given value '{val}', which is not an integer"
            )

        if defaults is not None:
            default_resources = cls.default(defaults)
            default_resources.update(cls(result))
            return default_resources

        return cls(result)

    @classmethod
    def from_mapping(cls, mapping: Mapping[str, ValidResource] | Self) -> Self:
        """Initialize from a mapping of resource to values.

        No validation is done on the keys. Values are converted into ``Resource``s and
        will be validated by that class.

        Raises
        ======
        ResourceDuplicationError
            if multiple variants of the same sized resource are provided (e.g. mem,
            mem_mb)
        ResourceTypeError:
            if a human readable resource with an incorrect type is given (e.g. a string
            to a suffixed resource)
        ResourceValidationError:
            if a given resource is not a ``str``, ``int``, ``float``, ``None``, or a
            callable.
        WorkflowError
            if a given resource is of the human-readable group but cannot be parsed
        """
        if isinstance(mapping, cls):
            return mapping

        return cls({key: Resource(key, val) for key, val in mapping.items()})

    def update(self, other: Resources | Mapping[str, ValidResource]):
        if isinstance(other, Resources):
            self._data.update(other._data)
        else:
            self._data.update(Resources.from_mapping(other))

    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def items(self):
        return self._data.items()

    def unwrapped_items(self) -> Iterator[Tuple[str, str | int | None]]:
        for key, val in self._data.items():
            yield key, val.value

    def copy(self):
        """Return deepcopy of the resource object."""
        return copy.deepcopy(self)

    def get(self, item: str):
        result = self._data.get(item)
        if result is None:
            return Resource("", None)
        return result

    def setdefault(self, ix: str, val: ValidResource | Resource | None):
        if ix in self._data:
            return self._data[ix]
        self[ix] = val
        return self._data[ix]

    def expand_items(
        self,
        *,
        constraints: Mapping[str, Resource | int | None],
        evaluate: Callable[[Resource], Resource] | None,
        skip: Optional[Collection[str]] = None,
        expand_sized: bool = True,
    ) -> Iterable[Tuple[str, str | int | TBDString | None]]:
        """Evaluate and constrain resources then convert into a simple resource mapping.

        SizedResources, including mem_mb and disk_mb, are optionally expanded into their
        human-readable and mebibyte versions for direct consumption in executors.

        Arguments
        =========
        constraints
            Mapping of resources to use as constraints.
        evaluate
            Callable that takes a single resource as argument. It should call that
            resource's ``evaluate()`` method with appropriate kwargs.
        skip
            Optional set of resource names that should not be expanded.
        expand_sized
            Enable/disable the expansion of sized resources (e.g. ``mem_mb``) into their
            unsuffixed (``mem``) and mebibyte (``mem_mib``) versions. Defaults to
            ``True``.
        """
        skip = set() if skip is None else skip
        for resource in self.values():
            if resource.name in skip:
                continue

            if evaluate is not None and resource.is_evaluable():
                resource = evaluate(resource)

            resource = resource.constrain(constraints.get(resource.name))

            yield resource.name, resource.value

            if (
                expand_sized
                and resource.name.endswith("_mb")
                and resource.name[:-3] in SizedResources
                and resource.value is not None
            ):
                yield resource.name[:-3], resource.format_human_friendly()
                mib = resource.from_mb_to_mib()
                yield mib.name, mib.value


ValidScope: TypeAlias = Literal["local"] | Literal["global"] | Literal["excluded"]


class ResourceScopes(Dict[str, ValidScope]):
    """Index of resource scopes, where each entry is 'RESOURCE': 'SCOPE'

    Each resource may be scoped as local, global, or excluded. Any resources not
    specified are considered global.
    """

    def __init__(
        self,
        *args: Mapping[str, ValidScope] | Tuple[str, ValidScope],
        **kwargs: ValidScope,
    ):
        super().__init__(*args, **kwargs)
        valid_scopes: set[ValidScope] = {"local", "global", "excluded"}
        for res in SizedResources:
            self._normalize_sizes(res)
        if set(self.values()) - valid_scopes:
            invalid_res = [
                res for res, scope in self.items() if scope not in valid_scopes
            ]
            invalid_pairs = {res: self[res] for res in invalid_res}

            # For now, we don't want excluded in the documentation
            raise ResourceScopesException(
                "Invalid resource scopes: entries must be defined as RESOURCE=SCOPE "
                "pairs, where SCOPE is either 'local' or 'global'",
                invalid_pairs,
            )

    def _normalize_sizes(self, resource: str):
        found: set[str] = set()
        for suffix in ["", "_mb", "_mib"]:
            if resource + suffix in self:
                found.add(resource + suffix)
        if not found or resource in found:
            return
        if len(found) > 1:
            raise ResourceDuplicationError(list(found))
        self[resource] = self.pop(found.pop())

    @classmethod
    def defaults(cls):
        return cls(mem_mb="local", disk_mb="local", runtime="excluded")

    def _test(self, resource: str, scope: ValidScope):
        for suffix in {"_mb", "_mib"}:
            stripped = resource.removesuffix(suffix)
            if stripped in SizedResources:
                resource = stripped
                break
        return self.get(resource) == scope

    def is_local(self, resource: str):
        """Resources are not tallied by the global scheduler when submitting jobs

        Each submitted job or group gets its own pool of the resource, as
        specified under --resources.


        Returns
        -------
        set
        """
        return self._test(resource, "local")

    def is_global(self, resource: str):
        """Resources tallied across all job and group submissions.

        Returns
        -------
        set
        """
        return self._test(resource, "global")

    def is_excluded(self, resource: str):
        """Resources not submitted to cluster jobs

        These resources are used exclusively by the global scheduler. The primary case
        is for additive resources in GroupJobs such as runtime, which would not be
        properly handled by the scheduler in the sub-snakemake instance. This scope is
        not currently intended for use by end-users and is thus not documented

        Returns
        -------
        set
        """
        return self._test(resource, "excluded")

    @property
    def locals(self):
        """Resources are not tallied by the global scheduler when submitting jobs

        Each submitted job or group gets its own pool of the resource, as
        specified under --resources.


        Returns
        -------
        set
        """
        return set(res for res, scope in self.items() if scope == "local")

    @property
    def globals(self):
        """Resources tallied across all job and group submissions.

        Returns
        -------
        set
        """
        return set(res for res, scope in self.items() if scope == "global")

    @property
    def excluded(self):
        """Resources not submitted to cluster jobs

        These resources are used exclusively by the global scheduler. The primary case
        is for additive resources in GroupJobs such as runtime, which would not be
        properly handled by the scheduler in the sub-snakemake instance. This scope is
        not currently intended for use by end-users and is thus not documented

        Returns
        -------
        set
        """
        suffixed_resources: set[str] = set()
        for resource in SizedResources:
            if resource in self:
                for suffix in {"_mb", "_mib"}:
                    suffixed_resources.add(resource + suffix)
        return (
            set(res for res, scope in self.items() if scope == "excluded")
            | suffixed_resources
        )

    def update(self, other: ResourceScopes | Mapping[str, ValidScope]):
        if not isinstance(other, ResourceScopes):
            other = ResourceScopes(other)
        super().update(other)


def is_ordinary_string(val: str):
    r"""Check if a string is an ordinary string.

    Ordinary strings are not evaluated and are not
    expected to be python expressions and be returned as is.

    Additionally, strings representing function calls, dictionary literals,
    or lambda expressions are considered offending strings.
    This function is useful for determining if a string can be safely
    returned to represent a resource value without further evaluation.
    It is important to note that this function does not check if the string
    is a valid Python identifier or a valid expression. It only checks
    if the string is an instance of `str` and does not match certain
    patterns that indicate it is a Python expression or callable.

    An ordinary string is defined as a string that:
    - Is an instance of the `str` type.
    - Does not match patterns that indicate it is a Python expression or a callable.

    The regular expression used for validation:
    - `^[a-zA-Z_]\w*\(.*\)$`: Matches function calls (e.g., `func_name(...)`).
    - `^\{.*\}$`: Matches strings that look like dictionary literals (e.g., `{...}`).
    - `^lambda\s.*:.*$`: Matches lambda expressions (e.g., `lambda x: x + 1`).
    - `.*[\+\-\*/\%].*|.*\.\w+.*`: Matches strings containing math operators or attribute access

    Parameters:
        val (any): The value to check.

    Returns:
        bool: True if the value is an ordinary string in this sense, False otherwise.
    """
    return isinstance(val, str) and not re.match(
        r"^[a-zA-Z_]\w*\(.*\)$|^\{.*\}$|^lambda\s.*:.*$|.*[\+\-\*/\%].*|.*\.\w+.*", val
    )
