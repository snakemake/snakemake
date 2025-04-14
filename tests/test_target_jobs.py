__authors__ = ["Tim Booth", "Johannes Köster"]
__copyright__ = "Copyright 2025, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from snakemake.target_jobs import parse_target_jobs_cli_args, TargetSpec

def test_parse_target_jobs():

    def ptj(args):
        return list(parse_target_jobs_cli_args(args))

    # Test empty list/None
    assert ptj(None) == []
    assert ptj([]) == []

    # Test a standard looking list of jobs in the old format.
    args1 = [ "ruleA:WILDCARD1=VALUE,WILDCARD2=VALUE",
              "ruleB:" ]
    assert ptj(args1) == \
        [ TargetSpec("ruleA", {"WILDCARD1": "VALUE",
                               "WILDCARD2": "VALUE"}),
          TargetSpec("ruleB", {}) ]

    # Test the version that uses JSON for robust quoting.
    args3 = [ 'ruleA::{ "WILDCARD1": "VALUE1", "WILDCARD2": "VALUE2" }',
              'ruleB::',
              'ruleA::{ "WILDCARD1": "VALUE3", "WILDCARD2": "VALUE4" }' ]
    assert ptj(args3) == \
        [ TargetSpec("ruleA", {"WILDCARD1": "VALUE1",
                               "WILDCARD2": "VALUE2"}),
          TargetSpec("ruleB", {}),
          TargetSpec("ruleA", {"WILDCARD1": "VALUE3",
                               "WILDCARD2": "VALUE4"}) ]

