from collections import defaultdict
import itertools as it
import operator as op
import re
import tempfile

from snakemake.exceptions import WorkflowError
from snakemake.common import TBDString


class DefaultResources:
    defaults = {
        "mem_mb": "max(2*input.size_mb, 1000)",
        "disk_mb": "max(2*input.size_mb, 1000)",
        "tmpdir": "system_tmpdir",
    }

    bare_defaults = {
        "tmpdir": "system_tmpdir",
    }

    @classmethod
    def decode_arg(cls, arg):
        try:
            return arg.split("=")
        except ValueError:
            raise ValueError("Resources have to be defined as name=value pairs.")

    @classmethod
    def encode_arg(cls, name, value):
        return "{}={}".format(name, value)

    def __init__(self, args=None, from_other=None, mode="full"):
        if mode == "full":
            self._args = dict(DefaultResources.defaults)
        elif mode == "bare":
            self._args = dict(DefaultResources.bare_defaults)
        else:
            raise ValueError("Unexpected mode for DefaultResources: {}".format(mode))

        if from_other is not None:
            self._args = dict(from_other._args)
            self.parsed = dict(from_other.parsed)
        else:
            if args is None:
                args = []

            self._args.update(
                {name: value for name, value in map(self.decode_arg, args)}
            )

            def fallback(val):
                def callable(wildcards, input, attempt, threads, rulename):
                    try:
                        value = eval(
                            val,
                            {
                                "input": input,
                                "attempt": attempt,
                                "threads": threads,
                                "system_tmpdir": tempfile.gettempdir(),
                            },
                        )
                    # Triggers for string arguments like n1-standard-4
                    except NameError:
                        return val
                    except Exception as e:
                        if not (
                            isinstance(e, FileNotFoundError) and e.filename in input
                        ):
                            # Missing input files are handled by the caller
                            raise WorkflowError(
                                "Failed to evaluate default resources value "
                                "'{}'.\n"
                                "    String arguments may need additional "
                                "quoting. Ex: --default-resources "
                                "\"tmpdir='/home/user/tmp'\".".format(val),
                                e,
                            )
                        raise e
                    return value

                return callable

            self.parsed = dict(_cores=1, _nodes=1)
            self.parsed.update(parse_resources(self._args, fallback=fallback))

    def set_resource(self, name, value):
        self._args[name] = "{}".format(value)
        self.parsed[name] = value

    @property
    def args(self):
        return [self.encode_arg(name, value) for name, value in self._args.items()]

    def __bool__(self):
        return bool(self.parsed)


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
                    # This will be caused by an inner job, which needs files created by
                    # the same group. All we can do is to ignore such jobs for now.
                    continue

            # Jobs in pipe groups must be run simultaneously, so we merge all the
            # resources of each pipe group into one big "job". Resources are combined
            # as "intralayer", so additives (like runtime) are maxed, and the rest are
            # summed
            for pipe in pipe_resources.values():
                job_resources.append(
                    cls._merge_resources(
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

            layers = cls._get_layers(int_resources, sorted_constraints.values(), sortby)

            # In each layer, sum or max across all resource types within the layer,
            # similar to summing along axis 0 in numpy:
            # [
            #   ( 3 ^ , 4 ^ , 1 ),
            #   ( 2 | , 1 | , 6 ),
            #   ( 1 | , 4 | , 0 ),
            # ]
            # The method is selected based on the label in int_resources and whether
            # that resource is additive
            intralayer_merge = [
                max if res in additive_resources else sum for res in int_resources
            ]
            combined = [
                [method(r) for method, r in zip(intralayer_merge, zip(*layer))]
                for layer in layers
            ]

            # Reassign the combined values from each layer to their resource names
            resorted = dict(zip(int_resources, zip(*combined)))

            for res, values in resorted.items():
                method = sum if res in additive_resources else max
                block_resources[res] = method(values)

            blocks.append(block_resources)

        if run_local:
            return {**cls._merge_resources(blocks, default_method=sum), "_nodes": 1}

        return {
            **cls._merge_resources(
                blocks,
                default_method=max,
                methods={res: sum for res in additive_resources},
            ),
            "_nodes": 1,
        }

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
    def _merge_resources(cls, resources, skip=[], methods={}, default_method=max):
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
                # TODO: This error is overly specific: if additive resources is ever
                # altered, the error will not remain correct.
                raise WorkflowError(
                    "Not enough resources were provided. This error is typically"
                    "caused by a Pipe group requiring too many resources. Note"
                    "that all resources including 'cores' are summed across"
                    "every member of the pipe group, except for runtime, which is"
                    "calculated via max().\n\n"
                    "Excess Resources:\n" + "\n".join(error_text)
                )

        # Remove final empty row. (The above loop ends each cycle by ensuring
        # there's an empty row)
        rows.pop()
        return rows


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
                    "Resource definition must start with a valid identifier, but found {}.".format(
                        res
                    )
                )
            try:
                val = int(val)
            except ValueError:
                if fallback is not None:
                    val = fallback(val)
                else:
                    raise ValueError(
                        "Resource definiton must contain an integer after the identifier."
                    )
            if res == "_cores":
                raise ValueError(
                    "Resource _cores is already defined internally. Use a different name."
                )
            resources[res] = val
    return resources
