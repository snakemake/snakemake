from collections import namedtuple
import typing

from snakemake_interface_executor_plugins.utils import TargetSpec

from snakemake.common import parse_key_value_arg


def parse_target_jobs_cli_args(target_jobs_args):
    errmsg = "Invalid target wildcards definition: entries have to be defined as WILDCARD=VALUE pairs"
    if target_jobs_args is not None:
        target_jobs = list()
        for entry in target_jobs_args:
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
