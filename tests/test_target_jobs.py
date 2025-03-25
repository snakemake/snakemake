__authors__ = ["Tim Booth", "Johannes Köster"]
__copyright__ = "Copyright 2025, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from snakemake.target_jobs import parse_target_jobs_cli_args, TargetSpec

def test_parse_target_jobs():

    # Test empty list/None
    assert parse_target_jobs_cli_args(None) is None
    assert parse_target_jobs_cli_args([]) == []

    # Test a standard looking list of jobs in the old format.
    args1 = [ "ruleA:WILDCARD1=VALUE,WILDCARD2=VALUE",
              "ruleB:" ]
    assert parse_target_jobs_cli_args(args1) == \
        [ TargetSpec("ruleA", {"WILDCARD1": "VALUE",
                               "WILDCARD2": "VALUE"}),
          TargetSpec("ruleB", {}) ]

    # Test parsing the same kind of thing in the new format
    args2 = [ "ruleA::WILDCARD1=VALUE,WILDCARD2=VALUE",
              "ruleB::" ]
    assert parse_target_jobs_cli_args(args2) == \
        [ TargetSpec("ruleA", {"WILDCARD1": "VALUE,WILDCARD2=VALUE"}),
          TargetSpec("ruleB", {}) ]

    # Test the version that allows for one arg per wildcard, to
    # make robust quoting easier
    args3 = [ "ruleA::WILDCARD1=VALUE",
              "ruleB:",
              "ruleB::",
              "ruleA::WILDCARD2=VALUE" ]
    assert parse_target_jobs_cli_args(args3) == \
        [ TargetSpec("ruleA", {"WILDCARD1": "VALUE",
                               "WILDCARD2": "VALUE"}),
          TargetSpec("ruleB", {}) ]

