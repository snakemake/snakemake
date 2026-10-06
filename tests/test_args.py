__authors__ = ["K.D. Murray"]
__copyright__ = "Copyright 2024, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


from snakemake.settings.types import Batch


def test_parse_batch():
    from snakemake.cli import parse_batch

    assert parse_batch("aggregate=1/2") == Batch("aggregate", 1, 2)


def test_target_jobs_wildcard_roundtrip():
    from snakemake.target_jobs import parse_target_jobs_cli_args
    from snakemake_interface_executor_plugins.utils import (
        TargetSpec,
        encode_target_jobs_cli_args,
    )

    # Quotes and commas must survive the
    # encode_target_jobs_cli_args -> parse_target_jobs_cli_args round trip.
    for want in [
        "'quoted name'",
        '"dq"',
        '5" pipe',
        "plain",
        "a,b",
        "'a,b'",
        "",
    ]:
        args = encode_target_jobs_cli_args([TargetSpec("r", {"w": want})])
        got = parse_target_jobs_cli_args(args)[0].wildcards_dict["w"]
        assert got == want, (want, args, got)

    args = encode_target_jobs_cli_args([TargetSpec("r", {"a": "1,2", "b": "'x'"})])
    assert parse_target_jobs_cli_args(args)[0].wildcards_dict == {
        "a": "1,2",
        "b": "'x'",
    }

    # A literal quote in an unwrapped (comma-free) value must not flip the
    # parser into quote mode and swallow the separator.
    args = encode_target_jobs_cli_args([TargetSpec("r", {"x": '5" pipe', "y": "next"})])
    assert parse_target_jobs_cli_args(args)[0].wildcards_dict == {
        "x": '5" pipe',
        "y": "next",
    }
