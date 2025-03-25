from collections import namedtuple, defaultdict
import typing

from snakemake_interface_executor_plugins.utils import TargetSpec

from snakemake.common import parse_key_value_arg

def parse_target_jobs_cli_args(target_jobs_args):
    """Parses a list of target job specs (--target-jobs) which are in the form
        RULE:WILDCARD1=VALUE,WILDCARD2=VALUE,...

        or

        RULE::WILDCARD1=VALUE RULE::WILDCARD2=VALUE

       Returns a list of TargetSpec objects.
    """
    errmsg = "Invalid target wildcards definition: entries have to be defined as WILDCARD=VALUE pairs"
    def parse_wildcard(entry, **kwargs):
        return parse_key_value_arg(entry, errmsg, **kwargs)

    if target_jobs_args is None:
        return None

    target_jobs = defaultdict(dict)
    for entry in target_jobs_args:
        rulename, wildcards = entry.split(":", 1)
        if wildcards == "" or wildcards == ":":
            target_jobs[rulename].update({})
        elif wildcards.startswith(":"):
            # New format indicated by double-colon.
            # Only one WILDCARD2=VALUE pair is allowed, and all quotes etc. will be kept.
            target_jobs[rulename].update([parse_wildcard(wildcards[1:], strip_quotes=False)])

        else:
            wc_dict = dict(
                parse_wildcard(entry) for entry in wildcards.split(",")
            )
            target_jobs[rulename].update(wc_dict)
    return [ TargetSpec(k, v) for k, v in target_jobs.items() ]
