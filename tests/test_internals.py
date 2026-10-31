import asyncio
import os
from collections import Counter
from types import SimpleNamespace

from snakemake.cli import parse_set_resources
from io import StringIO
from snakemake.profiles import ProfileConfigFileParser
import textwrap
import pytest

import snakemake.common.misc as common
from snakemake.common.misc import is_local_file
from snakemake.pathvars import Pathvars
from snakemake.ioutils import subpath
from snakemake.iocontainers import Wildcards, InputFiles
from snakemake.io import _IOFile
from snakemake.scheduling.job_scheduler import JobScheduler
from snakemake_interface_common.exceptions import WorkflowError


def test_subpath():
    assert subpath("test.txt", strip_suffix=".txt") == "test"
    assert subpath("test.txt", basename=True) == "test.txt"
    assert subpath("test.txt", parent=True) == "."
    assert subpath("test.txt", ancestor=1) == "."
    assert subpath("test.txt", ancestor=2) == "."

    assert (
        subpath(lambda wildcards, input: "test.txt", strip_suffix=".txt")(
            Wildcards(), input=InputFiles()
        )
        == "test"
    )
    assert subpath("results/foo/test.txt", parent=True) == "results/foo"
    assert subpath("results/foo/test.txt", ancestor=1) == "results/foo"
    assert subpath("results/foo/test.txt", ancestor=2) == "results"
    assert subpath("results/foo/test.txt", ancestor=3) == "."


def test_is_local_file_detects_local_paths():
    assert is_local_file("Snakefile")
    assert is_local_file("/tmp/Snakefile")
    assert is_local_file("file:///tmp/Snakefile")
    assert is_local_file("file:Snakefile")


def test_is_local_file_detects_windows_drive_paths(monkeypatch):
    monkeypatch.setattr(common, "ON_WINDOWS", True)

    assert is_local_file("C:/workflow/Snakefile")
    assert is_local_file("C:\\workflow\\Snakefile")


def test_is_local_file_rejects_remote_and_shorthand_uris():
    assert not is_local_file("https://example.com/workflow/Snakefile")
    assert not is_local_file("gh:owner/repo@main")
    assert not is_local_file("gl:group/project@main")
    assert not is_local_file("s3://bucket/workflow/Snakefile")


def test_pathvars_missing():
    with pytest.raises(WorkflowError):
        Pathvars.from_config({"pathvars": {"results": "myresults"}}).apply(
            "<bar>/log.txt"
        )


def test_pathvars_level():
    # test that pathvars from different levels have the correct priority
    # rule > module_config > module > config > workflow > default

    # default
    pv = Pathvars.with_defaults()
    assert pv.get("results") == "results"
    assert pv.get("logs") == "logs"
    assert pv.get("resources") == "resources"
    assert pv.get("benchmarks") == "benchmarks"
    assert pv.apply("<results>/foo.txt") == "results/foo.txt"
    assert pv.apply("<logs>/foo.txt") == "logs/foo.txt"
    assert pv.apply("<resources>/foo.txt") == "resources/foo.txt"
    assert pv.apply("<benchmarks>/foo.txt") == "benchmarks/foo.txt"

    # workflow
    pv.update(Pathvars.from_workflow({"results": "wf_results"}))
    assert pv.get("results") == "wf_results"
    assert pv.get("logs") == "logs"

    # config
    pv.update(Pathvars.from_config({"pathvars": {"logs": "cfg_logs"}}))
    assert pv.get("results") == "wf_results"
    assert pv.get("logs") == "cfg_logs"

    # module
    pv.update(Pathvars.from_module({"results": "mod_results", "logs": "mod_logs"}))
    assert pv.get("results") == "mod_results"
    assert pv.get("logs") == "mod_logs"

    # module_config
    pv.update(
        Pathvars.from_config(
            {"pathvars": {"logs": "modcfg_logs", "benchmarks": "modcfg_benchmarks"}},
            module_level=True,
        )
    )
    assert pv.get("results") == "mod_results"
    assert pv.get("logs") == "modcfg_logs"
    assert pv.get("benchmarks") == "modcfg_benchmarks"

    # rule
    pv.update(
        Pathvars.from_rule({"logs": "rule_logs", "benchmarks": "rule_benchmarks"})
    )
    assert pv.get("results") == "mod_results"
    assert pv.get("logs") == "rule_logs"
    assert pv.get("benchmarks") == "rule_benchmarks"
    assert pv.apply("<results>/foo.txt") == "mod_results/foo.txt"
    assert pv.apply("<logs>/foo.txt") == "rule_logs/foo.txt"


def test_profile_parse():
    profile = textwrap.dedent("""
    set-resources:
        rule1:
            slurm_partition: "42"
    """)
    stream = StringIO(profile)
    stream.name = "foo/config.yaml"

    parsed_profile = ProfileConfigFileParser().parse(stream)
    assert parsed_profile == {"set-resources": ["rule1:slurm_partition='42'"]}

    parsed_resources = parse_set_resources(parsed_profile["set-resources"])
    # ensure that resource is still a string after evaluation
    assert (
        parsed_resources["rule1"]["slurm_partition"]
        .evaluate(None, None, None, None, None, None)
        .value
        == "42"
    )


class _CountingIOFile(_IOFile):
    """_IOFile that counts the existence checks made while measuring its size."""

    calls = Counter()

    async def exists(self):
        _CountingIOFile.calls[str(self)] += 1
        return os.path.exists(self)

    async def size(self):
        return os.path.getsize(self)


def _update_input_sizes(scheduler, tmp_path, *job_inputs):
    # distinct _IOFile objects for the same path, as different jobs hold their own
    jobs = [
        SimpleNamespace(input=[_CountingIOFile(str(tmp_path / f), None) for f in files])
        for files in job_inputs
    ]
    asyncio.run(JobScheduler.update_input_sizes(scheduler, jobs))


def test_update_input_sizes_measures_existing_inputs_once(tmp_path):
    _CountingIOFile.calls = Counter()
    scheduler = SimpleNamespace(_input_sizes={})
    (tmp_path / "a").write_text("abc")
    (tmp_path / "empty").write_text("")

    _update_input_sizes(scheduler, tmp_path, ["a", "empty"], ["a", "missing"])
    (tmp_path / "missing").write_text("xy")
    _update_input_sizes(scheduler, tmp_path, ["a", "empty"], ["a", "missing"])

    # inputs that exist (even empty ones) are not measured again,
    # inputs that were missing are measured again
    assert _CountingIOFile.calls == {
        str(tmp_path / "a"): 1,
        str(tmp_path / "empty"): 1,
        str(tmp_path / "missing"): 2,
    }
    assert scheduler._input_sizes == {
        str(tmp_path / "a"): 3,
        str(tmp_path / "empty"): 0,
        str(tmp_path / "missing"): 2,
    }


def test_update_input_sizes_eviction_keeps_current_inputs(tmp_path):
    _CountingIOFile.calls = Counter()
    scheduler = SimpleNamespace(_input_sizes={f"old{i}": 1 for i in range(9999)})
    (tmp_path / "a").write_text("abc")
    (tmp_path / "b").write_text("xy")

    _update_input_sizes(scheduler, tmp_path, ["a"])
    # exceeds the cache limit: only sizes of the current inputs are kept,
    # including the one that is not measured again
    _update_input_sizes(scheduler, tmp_path, ["a"], ["b"])

    assert scheduler._input_sizes == {str(tmp_path / "a"): 3, str(tmp_path / "b"): 2}
    assert _CountingIOFile.calls == {str(tmp_path / "a"): 1, str(tmp_path / "b"): 1}
