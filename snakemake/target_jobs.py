from collections import defaultdict

from snakemake.common import parse_key_value_arg


def parse_target_jobs_cli_args(args):
    errmsg = "Invalid target wildcards definition: entries have to be defined as WILDCARD=VALUE pairs"
    if args.target_jobs is not None:
        target_jobs = defaultdict(dict)
        for entry in args.target_jobs:
            rulename, wildcards = entry.split(":", 1)
            if wildcards:
                for entry in wildcards.split(","):
                    wildcard, value = parse_key_value_arg(entry, errmsg=errmsg)
                    target_jobs[rulename][wildcard] = value
            else:
                target_jobs[rulename] = {}
        return target_jobs


def encode_target_jobs_cli_args(target_jobs):
    items = []
    for rulename, wildcards in target_jobs.items():
        wildcards = ",".join(f"{key}={value}" for key, value in wildcards.items())
        items.append(f"{rulename}:{wildcards}")
    return items
