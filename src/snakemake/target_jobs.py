from collections import namedtuple
import typing
import json

from snakemake_interface_executor_plugins.utils import TargetSpec

from snakemake.common import parse_key_value_arg

def parse_target_jobs_cli_args(target_jobs_args):
    """Parses a list of target job specs (--target-jobs) which are in the form
        RULE:WILDCARD1=VALUE,WILDCARD2=VALUE,...

        or, in the new format...

        RULE::<JSON-encoded dict>

       Returns a list of TargetSpec objects.
    """
    errmsg = "Invalid target wildcards definition: entries have to be defined as WILDCARD=VALUE pairs"
    def parse_wildcard(entry, **kwargs):
        return parse_key_value_arg(entry, errmsg, **kwargs)

    if target_jobs_args is None:
        return []

    # Note that we have to return a real list, not a generator
    target_list = []
    for entry in target_jobs_args:
        rulename, wildcards = entry.split(":", 1)
        if wildcards == "" or wildcards == ":":
            target_list.append(TargetSpec(rulename, {}))
        elif wildcards.startswith(":"):
            # New format indicated by double-colon. Expect JSON.
            target_list.append(TargetSpec(rulename, json.loads(wildcards[1:])))
        else:
            wc_dict = dict(
                parse_wildcard(entry) for entry in wildcards.split(",")
            )
            target_list.append(TargetSpec(rulename, wc_dict))

    return target_list
