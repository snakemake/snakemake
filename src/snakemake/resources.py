from __future__ import annotations
from collections import defaultdict
import copy
from dataclasses import dataclass
from humanfriendly import InvalidTimespan, InvalidSize, parse_size, parse_timespan
import itertools as it
import operator as op
import re
import shutil
import tempfile
import math
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Iterator,
    List,
    Literal,
    Mapping,
    Optional,
    Self,
    Sequence,
    Tuple,
    TypeAlias,
    cast,
)

from snakemake.common import get_input_function_aux_params, mb_to_mib, mib_to_mb
from snakemake.exceptions import (
    ResourceConstraintError,
    ResourceDuplicationError,
    ResourceError,
    ResourceInsufficiencyError,
    ResourceScopesException,
    ResourceValidationError,
    WorkflowError,
    is_file_not_found_error,
)
from snakemake.common.tbdstring import TBDString
from snakemake.io import AnnotatedString

if TYPE_CHECKING:
    from snakemake.jobs import Job
    from snakemake.io import Wildcards


@dataclass
class ParsedResource:
    orig_arg: str
    value: Any


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
                        for k, res, in job.resources.items()
                        if not isinstance(res, TBDString)
                    }
                    if job.pipe_group:
                        pipe_resources[job.pipe_group].append(res)
                    else:
                        job_resources.append(res)
                except FileNotFoundError:
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
                raise ResourceInsufficiencyError(additive_resources, err.args[0])

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
            return {**cls._merge_resource_dict(blocks, default_method=sum), "_nodes": 1}

        return {
            **cls._merge_resource_dict(
                blocks,
                default_method=max,
                methods={res: sum for res in additive_resources},
            ),
            "_nodes": 1,
        }

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
                    "Resource {name} is a string but not all group jobs require the "
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


def evaluable_from_mb_to_mib(name: str, func: Callable[..., Any]):
    def inner(*args: Any, **kwargs: Any):
        result = func(*args, **kwargs)
        if not isinstance(result, (int, float)):
            errmsg = (
                f"Evaluable resource must return a 'int' or 'float' to convert to mib. "
                f"{name} == {result} (type {type(result)})"
            )
            raise TypeError(errmsg)

        return mb_to_mib(round(result))

    return inner


def evaluable_from_mib_to_mb(name: str, func: Callable[..., Any]):
    def inner(*args: Any, **kwargs: Any):
        result = func(*args, **kwargs)
        if not isinstance(result, (int, float)):
            errmsg = (
                f"Evaluable resource must return a 'int' or 'float' to convert to mb. "
                f"{name} == {result} (type {type(result)})"
            )
            raise TypeError(errmsg)
        return mib_to_mb(round(result))

    return inner


ValidResource: TypeAlias = int | str | float | None | Callable[..., "ValidResource"]

SizedResources = {"mem", "disk"}
TimeResources = {"runtime"}
HumanFriendlyResources = TimeResources | SizedResources


class _NotComputedClass:
    pass


_NotComputed = _NotComputedClass()


class Resource:
    """Standardized representation of a resource name and value.

    In addition to ``int`` and ``str``, the standard resource types, ``Resource`` can
    store ``None``, signifying an unset resource, and a callable, signifying a resource
    to be calculated in the context of a rule. Via the ``.value`` property, ``Resource``
    ensures a resource has been computed before its value is retreived.

    The class additionally handles input validation and standardization.
    ``float``s get rounded into ``int``. If the resource name corresponds to a size
    resource ("disk" and "mem") or a timespan resource ("runtime"), strings values
    will be intepreted as a human readable resource and converted to an integer of
    appropriate denomination (Mb for size, hr for timespan).

    Arguments
    =========
    name: str
        Name of the resource (note that "disk", "mem", and "runtime" has special
        behaviour, as described above).
    value: int, str, Callable, None
        The concrete or unevaluated value
    """

    _evaluator: Callable[..., ValidResource] | None

    def __init__(self, name: str, value: ValidResource):
        if not (
            isinstance(value, (str, int, float)) or callable(value) or value is None
        ):
            raise ResourceValidationError(value)
        if isinstance(value, float):
            value = round(value)
        self.name = name
        try:
            self._value = self._parse_human_friendly(name, value)
        except InvalidSize as err:
            raise WorkflowError(
                f"Cannot parse '{name}' value into size in MB for setting "
                f"'{name}_mb' resource: {value}"
            ) from err
        except InvalidTimespan as err:
            raise WorkflowError(
                "Cannot parse 'runtime' value into minutes for setting 'runtime' "
                f"resource: {value}"
            ) from err
        if callable(self._value):
            self._evaluator = self._value
        elif isinstance(self._value, AnnotatedString) and self._value.is_callable():
            self._evaluator = cast(Callable[..., ValidResource], self._value.callable)
        else:
            self._evaluator = None

        if self._evaluator is None:
            self._evaluable = False
        else:
            self._evaluable = True

    def __repr__(self):
        if self.is_evaluable():
            value = "function(...)"
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
        return self._evaluable

    def evaluate(self, *args: Any, **kwargs: Any):
        """Evaluate the resource with given args.

        This method ensures the returned resource is correctly returned in a
        new ``Resource`` class, validated, and standardized. Resources should thus
        only be evaulated via this method.
        """
        if self._evaluator is None:
            return self
        kept_args = get_input_function_aux_params(self._evaluator, kwargs)
        return self.__class__(self.name, self._evaluator(*args, **kept_args))

    def constrain(self, other: Resource | int | None):
        """Use ``other`` as the maximum value for ``Self``, but only if both are integers.

        Returns self unaltered if either resource is ``None`` or if both are the same
        non-integer type.

        Errors if ``other`` has a different type from ``Self``, or if either are
        evaluable resources.
        """
        if isinstance(self._value, TBDString):
            return self

        other_val = other.value if isinstance(other, Resource) else other

        if self._value is None or other_val is None:
            return self

        self_val = self.value
        if not isinstance(self_val, type(other_val)):
            raise ResourceConstraintError(self_val, other_val)
        if isinstance(self_val, int):
            assert isinstance(other_val, int)
            return Resource(self.name, min(other_val, self_val))
        return self

    def to_mib(self):
        """Convert the resource to mebibytes.

        The name of the resource is updated to have a _mib suffix.

        The method does not track the current units of the resource, so will not prevent
        mebibytes from being "converted" into mebibytes multiple times.

        Un-evaluated resources can be converted.

        Errors if conversion is attempted on a str resource (or a callable that returns
        a str/None).
        """
        if isinstance(self._value, TBDString):
            value = self._value
        elif self.is_evaluable():
            value = evaluable_from_mb_to_mib(self.name, self._evaluator)
        elif not isinstance(self._value, int):
            errmsg = (
                f"Resource must be of type 'int' to convert to mib. {self.name} == "
                f"{self._value} (type {type(self._value)})"
            )
            raise TypeError(errmsg)
        else:
            value = mb_to_mib(self._value)
        return Resource(f"{self.name.removesuffix('_mb')}_mib", value)

    def to_mb(self):
        """Convert the resource to megabytes.

        The name of the resource is updated to have a _mb suffix.

        The method does not track the current units of the resource, so will not prevent
        mebibytes from being "converted" into megabytes multiple times.

        Un-evaluated resources can be converted.

        Errors if conversion is attempted on a str resource (or a callable that returns
        a str/None).
        """
        if isinstance(self._value, TBDString):
            value = self._value
        elif self.is_evaluable():
            value = evaluable_from_mib_to_mb(self.name, self._evaluator)
        elif not isinstance(self._value, int):
            errmsg = (
                f"Resource must be of type 'int' to convert to mb. {self.name} == "
                f"{self._value} (type {type(self._value)})"
            )
            raise TypeError(errmsg)
        else:
            value = mib_to_mb(self._value)
        return Resource(f"{self.name.removesuffix('_mib')}_mb", value)

    def without_suffix(self):
        """Converts a suffixed resource (e.g. mem_mb) into its unsuffixed version (mb).

        If the resource has the _mib suffix, the value is converted into megabytes.
        """
        if self.name.endswith("_mb"):
            return self.with_name(self.name.removesuffix("_mb"))
        if self.name.endswith("_mib"):
            return self.to_mb().with_name(self.name.removesuffix("_mb"))
        return self

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
            python expression to be evaulated
        with_threads_arg: boolean
            If True, include ``threads`` as an argument in the returned function
        """

        def threads_evaluator(
            wildcards: Wildcards,
            input: Any,
            attempt: int,
            threads: int,
            rulename: str,
        ):
            return cls.cli_evaluator(
                name,
                value,
                input=input,
                attempt=attempt,
                threads=threads,
            )

        def nonthreads_evaluator(
            wildcards: Wildcards, input: Any, attempt: int, rulename: str
        ):
            return cls.cli_evaluator(
                name,
                value,
                input=input,
                attempt=attempt,
            )

        if with_threads_arg:
            return Resource(name, threads_evaluator)

        return Resource(name, nonthreads_evaluator)

    @staticmethod
    def cli_evaluator(
        name: str, val: str, *, input: Any, attempt: int, threads: int | None = None
    ):
        """Evaluate a python expression to determine resource value.

        Intended to be used from ``from_cli_expression``.
        """
        try:
            value = eval(
                val,
                {
                    "input": input,
                    "attempt": attempt,
                    "system_tmpdir": tempfile.gettempdir(),
                    "shutil": shutil,
                    **({"threads": threads} if threads is not None else {}),
                },
            )
        # Triggers for string arguments like n1-standard-4
        except (NameError, SyntaxError):
            return val
        except Exception as e:
            try:
                # check if resource is parsable as human friendly (given the correct
                # name and formatted value). If it is, we return the parsed value to
                # save a step later.
                if name in HumanFriendlyResources:
                    return Resource._parse_human_friendly(name, val)
            except (InvalidSize, InvalidTimespan):
                pass
            if is_file_not_found_error(e, input):
                return TBDString()

            raise WorkflowError(
                "Failed to evaluate resources value "
                f"'{val}'.\n"
                "    String arguments may need additional "
                "quoting. E.g.: --default-resources "
                "\"tmpdir='/home/user/tmp'\" or "
                "--set-resources \"somerule:someresource='--nice=100'\". "
                "This also holds for setting resources inside of a profile, where "
                "you might have to enclose them in single and double quotes, "
                "i.e. someresource: \"'--nice=100'\".",
                e,
            )
        return value

    @staticmethod
    def _parse_human_friendly(name: str, value: ValidResource):
        if isinstance(value, str) and not isinstance(value, TBDString):
            value = value.strip("'\"")
            if name in SizedResources:
                return max(int(math.ceil(parse_size(value) / 1e6)), 1)
            elif name in TimeResources:
                return max(int(round(parse_timespan(value) / 60)), 1)
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
    """

    DEFAULTS = {
        "full": {
            "mem": "min(max(2*input.size_mb, 1000), 8000)",
            "disk": "max(2*input.size_mb, 1000)",
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
        for suffix in ["", "_mb", "_mib"]:
            if resource + suffix in self._data:
                found.add(resource + suffix)
        if not found or {resource} == found:
            return
        if len(found) > 1:
            raise ResourceDuplicationError(list(found))
        self._data[resource] = self._data.pop(found.pop()).without_suffix()

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
            f"{name}={value.value}"
            for name, value in self._data.items()
            if not value.is_evaluable() and value.value is not None
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
        only_positive_integers: bool = False,
        defaults: None | Literal["bare"] | Literal["full"] = None,
    ) -> Callable[[List[str]], Self]:
        """Return a parsing function with preset keyword arguments.

        Intended for use with argparse.
        """
        return lambda exprs: cls.parse(
            exprs,
            allow_expressions=allow_expressions,
            only_positive_integers=only_positive_integers,
            defaults=defaults,
        )

    @classmethod
    def parse(
        cls,
        exprs: List[str],
        *,
        allow_expressions: bool = False,
        only_positive_integers: bool = False,
        defaults: None | Literal["bare"] | Literal["full"] = None,
    ):
        """Parse a series of CLI-style string resource assigments.

        Arguments
        =========
        exprs: list of str
            Each expr is in the form RESOURCE=VALUE, where RESOURCE is a valid python
            identifier, and VALUE is a number, a quoted string, or a valid python
            expression.
        allow_expressions:
            Allows use of quoted strings and python expressions in the list of exprs.
            Otherwise, an error is raised if these are found.
        only_positive_integers:
            Requires integers in the list of exprs to be positive.
        defaults: "bare" or "full"
            Initializes the returned resource object with the given default set, if
            provided. Provided exprs will override default resources.
        """
        valid = re.compile(r"[a-zA-Z_]\w*$")
        unparsed = dict([cls.decode_arg(expr) for expr in exprs])

        result: Dict[str, Resource] = {}
        for res, val in unparsed.items():
            if res == "_cores":
                raise ValueError(
                    "Resource _cores is already defined internally. Use a different "
                    "name."
                )
            if not valid.match(res):
                raise ValueError(
                    "Resource definition must start with a valid identifier, but found "
                    "{}.".format(res)
                )

            try:
                val = int(val)
            except ValueError:
                pass
            else:
                if only_positive_integers and val < 0:
                    raise ValueError("integer resource must be positive")
                result[res] = Resource(res, val)
                continue

            if allow_expressions:
                result[res] = Resource.from_cli_expression(
                    res, val, with_threads_arg=True
                )
                continue

            raise ValueError(
                "Resource definition must contain an integer after the identifier"
            )

        if defaults is not None:
            default_resources = cls.default(defaults)
            default_resources.update(cls(result))
            return default_resources

        return cls(result)

    @classmethod
    def from_mapping(cls, mapping: Mapping[str, Any]):
        """Initialize from a mapping of resource to values.

        No validation is done on the keys. Values are converted into ``Resource``s and
        will be validated by that class.
        """
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

    def unwrapped_nonstr_items(self) -> Iterator[Tuple[str, int | None]]:
        for key, val in self._data.items():
            val = val.value
            if isinstance(val, str):
                continue
            yield key, val

    def copy(self):
        """Return deepcopy of the resource object."""
        return copy.deepcopy(self)

    def get(self, item: str):
        result = self._data.get(item)
        if result is None:
            return Resource("", None)
        return result


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
