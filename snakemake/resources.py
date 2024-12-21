from collections import UserDict, defaultdict
from dataclasses import dataclass
import itertools as it
import operator as op
import re
import tempfile
import math
from typing import Any

from snakemake.exceptions import (
    ResourceScopesException,
    WorkflowError,
    is_file_not_found_error,
)
from snakemake.common.tbdstring import TBDString
from snakemake.logging import logger


@dataclass
class ParsedResource:
    orig_arg: str
    value: Any


class DefaultResources:
    defaults = {
        "mem_mb": "min(max(2*input.size_mb, 1000), 8000)",
        "disk_mb": "max(2*input.size_mb, 1000)",
        "tmpdir": "system_tmpdir",
    }

    bare_defaults = {"tmpdir": "system_tmpdir"}

    @classmethod
    def decode_arg(cls, arg):
        try:
            return arg.split("=", maxsplit=1)
        except ValueError:
            raise ValueError("Resources have to be defined as name=value pairs.")

    @classmethod
    def encode_arg(cls, name, value):
        return f"{name}={value}"

    def __init__(self, args=None, from_other=None, mode="full"):
        if mode == "full":
            self._args = dict(DefaultResources.defaults)
        elif mode == "bare":
            self._args = dict(DefaultResources.bare_defaults)
        else:
            raise ValueError(f"Unexpected mode for DefaultResources: {mode}")

        if from_other is not None:
            self._args = dict(from_other._args)
            self.parsed = dict(from_other.parsed)
        else:
            if args is None:
                args = []

            self._args.update(
                {name: value for name, value in map(self.decode_arg, args)}
            )

            self.parsed = dict(_cores=1, _nodes=1)
            self.parsed.update(
                parse_resources(self._args, fallback=eval_resource_expression)
            )

    def set_resource(self, name, value):
        self._args[name] = f"{value}"
        self.parsed[name] = value

    @property
    def args(self):
        return [self.encode_arg(name, value) for name, value in self._args.items()]

    def __bool__(self):
        return bool(self.parsed)


class ResourceScopes(UserDict):
    """Index of resource scopes, where each entry is 'RESOURCE': 'SCOPE'

    Each resource may be scoped as local, global, or excluded. Any resources not
    specified are considered global.
    """

    def __init__(self, *args, **kwargs):
        self.data = dict(*args, **kwargs)
        valid_scopes = {"local", "global", "excluded"}
        if set(self.data.values()) - valid_scopes:
            invalid_res = [
                res for res, scope in self.data.items() if scope not in valid_scopes
            ]
            invalid_pairs = {res: self.data[res] for res in invalid_res}

            # For now, we don't want excluded in the documentation
            raise ResourceScopesException(
                "Invalid resource scopes: entries must be defined as RESOURCE=SCOPE "
                "pairs, where SCOPE is either 'local' or 'global'",
                invalid_pairs,
            )

    @classmethod
    def defaults(cls):
        return cls(mem_mb="local", disk_mb="local", runtime="excluded")

    @property
    def locals(self):
        """Resources are not tallied by the global scheduler when submitting jobs

        Each submitted job or group gets its own pool of the resource, as
        specified under --resources.


        Returns
        -------
        set
        """
        return set(res for res, scope in self.data.items() if scope == "local")

    @property
    def globals(self):
        """Resources tallied across all job and group submissions.

        Returns
        -------
        set
        """
        return set(res for res, scope in self.data.items() if scope == "global")

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
        return set(res for res, scope in self.data.items() if scope == "excluded")


class GroupResources:
    @classmethod
    def basic_layered(
        cls,
        toposorted_jobs,
        constraints,
        run_local,
        additive_resources=None,
        sortby=None,
    ):
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

        total_resources = defaultdict(int)
        total_resources["_nodes"] = 1
        blocks = []
        # iterate over siblings that can be executed in parallel
        for siblings in toposorted_jobs:
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
            int_resources = {}
            # Sort all integer resources in job_resources into int_resources. Resources
            # defined as a string are placed immediately into block_resources.
            for res in resource_types:
                if res == "_nodes":
                    continue

                values = [resources.get(res, 0) for resources in job_resources]

                if cls._is_string_resource(res, values):
                    block_resources[res] = values[0]
                else:
                    int_resources[res] = values

            # Collect values from global_resources to use as constraints.
            sorted_constraints = {
                name: constraints.get(name, None) for name in int_resources
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
                    int_resources, sorted_constraints.values(), sortby
                )
            except WorkflowError as err:
                raise cls._get_saturated_resource_error(additive_resources, err.args[0])

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
    def _get_saturated_resource_error(cls, additive_resources, excess_resources):
        isare = "is" if len(additive_resources) == 1 else "are"
        additive_clause = (
            (f", except for {additive_resources}, which {isare} calculated via max(). ")
            if additive_resources
            else ". "
        )
        return WorkflowError(
            "Not enough resources were provided. This error is typically "
            "caused by a Pipe group requiring too many resources. Note "
            "that resources are summed across every member of the pipe "
            f"group{additive_clause}"
            f"Excess Resources:\n{excess_resources}"
        )

    @classmethod
    def _is_string_resource(cls, name, values):
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
    def _merge_resource_dict(cls, resources, skip=[], methods={}, default_method=max):
        grouped = {}
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

        ret = {}
        for res, values in grouped.items():
            if res in skip:
                continue

            if cls._is_string_resource(res, values):
                ret[res] = values[0]
            elif res in methods:
                ret[res] = methods[res](values)
            else:
                ret[res] = default_method(values)
        return ret

    @classmethod
    def _merge_resource_layer(cls, resources, methods):
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
    def _check_constraint(resources, constraints):
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
    def _get_layers(cls, resources, constraints, sortby=None):
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

        # By zipping, we combine the vals into tuples based on job, 1 tuple per
        # job: [ (val1, 1_val1, 2_val1), ...]. In each tuple, the resources
        # will remain in the same order as the original dict, so their identity
        # can be extracted later.
        resource_groups = zip(*resources.values())

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
                raise WorkflowError("\n".join(error_text))

        # Remove final empty row. (The above loop ends each cycle by ensuring
        # there's an empty row)
        rows.pop()
        return rows


def eval_resource_expression(val, threads_arg=True):
    def generic_callable(val, threads_arg, **kwargs):
        args = {
            "input": kwargs["input"],
            "attempt": kwargs["attempt"],
            "system_tmpdir": tempfile.gettempdir(),
        }
        if threads_arg:
            args["threads"] = kwargs["threads"]
        try:
            value = eval(
                val,
                args,
            )
        # Triggers for string arguments like n1-standard-4
        except NameError:
            return val
        except Exception as e:
            if is_humanfriendly_resource(val):
                return val
            if not is_file_not_found_error(e, kwargs["input"]):
                # Missing input files are handled by the caller
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
            raise e
        return value

    if threads_arg:

        def callable(wildcards, input, attempt, threads, rulename):
            return generic_callable(
                val,
                threads_arg=threads_arg,
                wildcards=wildcards,
                input=input,
                attempt=attempt,
                threads=threads,
                rulename=rulename,
            )

    else:

        def callable(wildcards, input, attempt, rulename):
            return generic_callable(
                val,
                threads_arg=threads_arg,
                wildcards=wildcards,
                input=input,
                attempt=attempt,
                rulename=rulename,
            )

    return callable


def parse_resources(resources_args, fallback=None):
    """Parse resources from args."""
    resources = dict()

    if resources_args is not None:
        valid = re.compile(r"[a-zA-Z_]\w*$")

        if isinstance(resources_args, list):
            resources_args = map(DefaultResources.decode_arg, resources_args)
        else:
            resources_args = resources_args.items()

        for res, val in resources_args:
            if not valid.match(res):
                raise ValueError(
                    "Resource definition must start with a valid identifier, but found "
                    "{}.".format(res)
                )

            try:
                val = int(val)
            except ValueError:
                if fallback is not None:
                    val = fallback(val)
                else:
                    raise ValueError(
                        "Resource definition must contain an integer, string or python expression after the identifier."
                    )
            if res == "_cores":
                raise ValueError(
                    "Resource _cores is already defined internally. Use a different "
                    "name."
                )
            resources[res] = val
    return resources


def infer_resources(name, value, resources: dict):
    """Infer resources from a given one, if possible."""
    from humanfriendly import parse_size, parse_timespan, InvalidTimespan, InvalidSize

    if isinstance(value, str):
        value = value.strip("'\"")

    if (
        (name == "mem" or name == "disk")
        and isinstance(value, str)
        and not isinstance(value, TBDString)
    ):
        inferred_name = f"{name}_mb"
        try:
            in_bytes = parse_size(value)
        except InvalidSize:
            raise WorkflowError(
                f"Cannot parse mem or disk value into size in MB for setting {inferred_name} resource: {value}"
            )
        resources[inferred_name] = max(int(math.ceil(in_bytes / 1e6)), 1)
    elif (
        name == "runtime"
        and isinstance(value, str)
        and not isinstance(value, TBDString)
    ):
        try:
            parsed = max(int(round(parse_timespan(value) / 60)), 1)
        except InvalidTimespan:
            raise WorkflowError(
                f"Cannot parse runtime value into minutes for setting runtime resource: {value}"
            )
        resources["runtime"] = parsed


def is_humanfriendly_resource(value):
    from humanfriendly import parse_size, parse_timespan, InvalidTimespan, InvalidSize

    try:
        parse_size(value)
        return True
    except InvalidSize:
        pass

    try:
        parse_timespan(value)
        return True
    except InvalidTimespan:
        pass

    return False
