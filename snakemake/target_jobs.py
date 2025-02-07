from collections import namedtuple
import typing
import re
from snakemake_interface_executor_plugins.utils import TargetSpec

from snakemake.common import parse_key_value_arg


def parse_target_jobs_cli_args(target_jobs_args):
    errmsg = "Invalid target wildcards definition: entries have to be defined as WILDCARD=VALUE pairs"
    if target_jobs_args is not None:
        target_jobs = list()
        for entry in target_jobs_args:
            rulename, wildcards = entry.split(":", 1)
            if wildcards:

                def remove_len_zero(tpl):
                    return tuple(val for val in tpl if len(val)>0)

                matches = re.findall(r'(\w+)="([^"]*)",|(\w+)=([^,]*),', wildcards+",")
                wildcards = dict(map(remove_len_zero, matches))
                target_jobs.append(TargetSpec(rulename, wildcards))
                target_jobs.append(TargetSpec(rulename, wildcards))
            else:
                target_jobs.append(TargetSpec(rulename, dict()))
        return target_jobs
