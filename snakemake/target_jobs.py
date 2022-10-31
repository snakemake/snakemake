from collections import namedtuple
import typing

from snakemake.common import parse_key_value_arg


TargetSpec = namedtuple("TargetSpec", ["rulename", "wildcards_dict"])


def parse_target_jobs_cli_args(args):
    errmsg = "Invalid target wildcards definition: entries have to be defined as WILDCARD=VALUE pairs"
    if args.target_jobs is not None:
        target_jobs = list()
        for entry in args.target_jobs:
            rulename, wildcards = entry.split(":", 1)
            if wildcards:

                def parse_wildcard(entry):
                    return parse_key_value_arg(entry, errmsg)

                wildcards = dict(
                    parse_wildcard(entry) for entry in wildcards.split(",")
                )
                target_jobs.append(TargetSpec(rulename, wildcards))
            else:
                target_jobs.append(TargetSpec(rulename, dict()))
        return target_jobs


def encode_target_jobs_cli_args(
    target_jobs: typing.List[TargetSpec],
) -> typing.List[str]:
    items = []
    for spec in target_jobs:
        wildcards = ",".join(
            f"{key}={value}" for key, value in spec.wildcards_dict.items()
        )
        items.append(f"{spec.rulename}:{wildcards}")
    return items
