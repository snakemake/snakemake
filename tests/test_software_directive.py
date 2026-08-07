__authors__ = ["Yun Jiang"]
__copyright__ = "Copyright 2026, Yun Jiang"
__license__ = "MIT"

import os
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

sys.path.insert(0, os.path.dirname(__file__))

from snakemake_interface_software_deployment_plugins import EnvSpecSourceFile
from snakemake_software_deployment_plugin_conda import EnvSpec as CondaEnvSpec
from snakemake_software_deployment_plugin_container import EnvSpec as ContainerEnvSpec
from snakemake_software_deployment_plugin_envmodules import EnvSpec as EnvModuleEnvSpec

from snakemake.deployment import EnvSpecs, SoftwareDeploymentManager
from snakemake.sourcecache import _replace_suffix


def _make_conda_spec(path: str = "env.yaml") -> CondaEnvSpec:
    spec = CondaEnvSpec(envfile=EnvSpecSourceFile(path))
    spec.technical_init()
    return spec


def _make_container_spec(image: str = "ubuntu:22.04") -> ContainerEnvSpec:
    spec = ContainerEnvSpec(image)
    spec.technical_init()
    return spec


def _make_envmodules_spec(*names: str) -> EnvModuleEnvSpec:
    spec = EnvModuleEnvSpec(*names)
    spec.technical_init()
    return spec


# ---------------------------------------------------------------------------
# EnvSpecs.interpret() unit tests
# ---------------------------------------------------------------------------


class TestEnvSpecsInterpret:
    """Unit tests for EnvSpecs.interpret() — the core logic that converts
    legacy and new directives into a unified EnvSpecBase tree."""

    def test_empty_specs_returns_none(self):
        specs = EnvSpecs()
        # interpret() has an assert not is_callable(), then checks software_spec
        # When everything is None, it returns None
        result = specs.interpret()
        assert result is None

    def test_software_spec_priority_over_legacy(self):
        """When both software_spec and legacy_conda_env are set,
        software_spec takes priority."""
        new_spec = _make_conda_spec("new-env.yaml")
        specs = EnvSpecs(
            software_spec=new_spec,
            legacy_conda_env="legacy-env.yaml",
        )
        result = specs.interpret()
        assert result is new_spec

    def test_software_spec_only(self):
        spec = _make_conda_spec("env.yaml")
        specs = EnvSpecs(software_spec=spec)
        result = specs.interpret()
        assert result is spec
        assert result.kind == "conda"

    def test_legacy_conda_yaml_creates_conda_spec(self):
        specs = EnvSpecs(legacy_conda_env="env.yaml")
        result = specs.interpret()
        assert isinstance(result, CondaEnvSpec)
        assert result.kind == "conda"

    def test_legacy_conda_yml_creates_conda_spec(self):
        specs = EnvSpecs(legacy_conda_env="env.yml")
        result = specs.interpret()
        assert isinstance(result, CondaEnvSpec)

    def test_legacy_conda_name_creates_conda_spec(self):
        specs = EnvSpecs(legacy_conda_env="my-env-name")
        result = specs.interpret()
        assert isinstance(result, CondaEnvSpec)

    def test_legacy_container_creates_container_spec(self):
        specs = EnvSpecs(legacy_container_img="docker://ubuntu:22.04")
        result = specs.interpret()
        assert isinstance(result, ContainerEnvSpec)
        assert result.kind == "container"

    def test_legacy_envmodules_creates_envmodules_spec(self):
        specs = EnvSpecs(legacy_env_modules=["samtools", "bwa"])
        result = specs.interpret()
        assert isinstance(result, EnvModuleEnvSpec)
        assert result.kind == "envmodules"

    def test_conda_within_container(self):
        """Legacy pattern: both conda and container → conda within container."""
        specs = EnvSpecs(
            legacy_conda_env="env.yaml",
            legacy_container_img="ubuntu:22.04",
        )
        result = specs.interpret()
        # When both legacy_conda and legacy_container are set,
        # conda spec should be the primary, with within=container
        assert isinstance(result, CondaEnvSpec)
        assert result.within is not None
        assert isinstance(result.within, ContainerEnvSpec)

    def test_envmodules_fallback_to_conda(self):
        """Legacy pattern: envmodules + conda → envmodules with fallback=conda."""
        specs = EnvSpecs(
            legacy_env_modules=["samtools"],
            legacy_conda_env="env.yaml",
        )
        result = specs.interpret()
        assert isinstance(result, EnvModuleEnvSpec)
        assert result.fallback is not None
        assert isinstance(result.fallback, CondaEnvSpec)

    def test_envmodules_fallback_to_container(self):
        """Legacy pattern: envmodules + container → envmodules with fallback=container."""
        specs = EnvSpecs(
            legacy_env_modules=["samtools"],
            legacy_container_img="ubuntu:22.04",
        )
        result = specs.interpret()
        assert isinstance(result, EnvModuleEnvSpec)
        assert result.fallback is not None
        assert isinstance(result.fallback, ContainerEnvSpec)

    def test_conda_within_container_and_envmodules_fallback(self):
        """Triple legacy: envmodules + conda + container.
        envmodules is primary, fallback to conda-within-container."""
        specs = EnvSpecs(
            legacy_env_modules=["samtools"],
            legacy_conda_env="env.yaml",
            legacy_container_img="ubuntu:22.04",
        )
        result = specs.interpret()
        assert isinstance(result, EnvModuleEnvSpec)
        # fallback should be conda (which has within=container)
        assert isinstance(result.fallback, CondaEnvSpec)
        assert isinstance(result.fallback.within, ContainerEnvSpec)


class TestEnvSpecsProperties:
    """Tests for EnvSpecs helper methods."""

    def test_is_empty_all_none(self):
        specs = EnvSpecs()
        assert specs.is_empty()

    def test_is_empty_with_software_spec(self):
        spec = _make_conda_spec()
        specs = EnvSpecs(software_spec=spec)
        assert not specs.is_empty()

    def test_is_empty_with_legacy_conda(self):
        specs = EnvSpecs(legacy_conda_env="env.yaml")
        assert not specs.is_empty()

    def test_is_callable_with_callable_software_spec(self):
        specs = EnvSpecs(software_spec=lambda wc: "env.yaml")
        assert specs.is_callable()

    def test_is_callable_with_non_callable(self):
        spec = _make_conda_spec()
        specs = EnvSpecs(software_spec=spec)
        assert not specs.is_callable()

    def test_contains_conda_or_container_software_spec_conda(self):
        spec = _make_conda_spec()
        specs = EnvSpecs(software_spec=spec)
        assert specs.contains_conda_or_container()

    def test_contains_conda_or_container_software_spec_container(self):
        spec = _make_container_spec()
        specs = EnvSpecs(software_spec=spec)
        assert specs.contains_conda_or_container()

    def test_contains_conda_or_container_envmodules_only(self):
        spec = _make_envmodules_spec("samtools")
        specs = EnvSpecs(software_spec=spec)
        # envmodules alone does not contain conda or container
        assert not specs.contains_conda_or_container()

    def test_contains_conda_or_container_legacy_conda(self):
        specs = EnvSpecs(legacy_conda_env="env.yaml")
        assert specs.contains_conda_or_container()

    def test_contains_conda_or_container_legacy_container(self):
        specs = EnvSpecs(legacy_container_img="ubuntu:22.04")
        assert specs.contains_conda_or_container()

    def test_resolve_callables(self):
        specs = EnvSpecs(software_spec=lambda wc: "env.yaml")
        assert specs.is_callable()
        resolved = specs.resolve_callables(lambda fn: fn(None))
        assert not resolved.is_callable()


class TestEnvSpecComposition:
    """Tests for within and fallback composition on EnvSpecBase."""

    def test_within_set_correctly(self):
        conda_spec = _make_conda_spec()
        container_spec = _make_container_spec()
        conda_spec.within = container_spec
        assert conda_spec.within is container_spec
        assert conda_spec.within.kind == "container"

    def test_fallback_set_correctly(self):
        envmodules_spec = _make_envmodules_spec("samtools")
        conda_spec = _make_conda_spec()
        envmodules_spec.fallback = conda_spec
        assert envmodules_spec.fallback is conda_spec
        assert envmodules_spec.fallback.kind == "conda"

    def test_nested_within_and_fallback(self):
        """envmodules → fallback=conda(within=container)."""
        conda_spec = _make_conda_spec()
        container_spec = _make_container_spec()
        conda_spec.within = container_spec

        envmodules_spec = _make_envmodules_spec("samtools")
        envmodules_spec.fallback = conda_spec

        assert envmodules_spec.fallback.within is container_spec
        assert envmodules_spec.fallback.kind == "conda"
        assert envmodules_spec.fallback.within.kind == "container"

    def test_has_source_paths_with_conda_envfile(self):
        spec = _make_conda_spec("env.yaml")
        assert spec.has_source_paths()

    def test_has_source_paths_container_only(self):
        spec = _make_container_spec()
        assert not spec.has_source_paths()

    def test_has_source_paths_nested_within(self):
        """Source path in nested within spec should propagate."""
        container_spec = _make_container_spec()
        conda_spec = _make_conda_spec("env.yaml")
        conda_spec.within = container_spec
        assert conda_spec.has_source_paths()


def _make_deployment_manager(*selected_methods: str) -> SoftwareDeploymentManager:
    """Build a SoftwareDeploymentManager with a stub workflow.

    The constructor only reads ``deployment_settings.deployment_methods``; the
    real plugin registry auto-discovers the installed conda/container/envmodules
    plugins. This lets us unit-test the ``selected_plugin_kinds`` gating in
    ``get_env`` (deployment/__init__.py) without spinning up a full workflow or
    performing any real deployment.
    """
    stub_workflow = SimpleNamespace(
        deployment_settings=SimpleNamespace(deployment_methods=set(selected_methods))
    )
    return SoftwareDeploymentManager(stub_workflow)


class TestSuffixReplacement:
    """Tests for the suffix-replacement mechanism (interface PR #54) that backs
    ``within=`` composition: a source file's suffix is swapped (e.g. to derive a
    manifest from a template). Covers the ``_replace_suffix`` helper that
    ``rules.py`` relies on, plus the ``EnvSpecSourceFile.suffix_replacement``
    field contract."""

    def test_replace_suffix_yaml(self):
        assert _replace_suffix("env.yaml", [".yaml", ".yml"], ".lock") == "env.lock"

    def test_replace_suffix_yml(self):
        assert _replace_suffix("env.yml", [".yaml", ".yml"], ".lock") == "env.lock"

    def test_replace_suffix_no_match_returns_none(self):
        assert _replace_suffix("env.txt", [".yaml", ".yml"], ".lock") is None

    def test_replace_suffix_multi_part_suffix(self):
        # multi-part suffixes like .tar.gz must match as a whole token
        assert _replace_suffix("pkg.tar.gz", [".tar.gz", ".zip"], ".tar.xz") == (
            "pkg.tar.xz"
        )

    def test_replace_suffix_first_match_wins(self):
        # candidate suffixes are tried in list order
        assert _replace_suffix("env.yaml", [".yml", ".yaml"], ".lock") == "env.lock"

    def test_envspecsourcefile_exposes_suffix_replacement_field(self):
        # the interface dataclass carries a suffix_replacement slot (PR #54);
        # it defaults to None so specs that do not use it are unaffected
        source_file = EnvSpecSourceFile("env.yaml")
        assert hasattr(source_file, "suffix_replacement")
        assert source_file.suffix_replacement is None


class TestCrossKindFallback:
    """Tests for ``SoftwareDeploymentManager.get_env`` cross-kind fallback, which
    gates env resolution on ``selected_plugin_kinds`` (the redesign in
    deployment/__init__.py). When a spec's kind is not among the activated
    plugins, resolution falls through to the spec's fallback (if any); this is
    what makes ``container(...) or conda(...)`` work when only one kind is
    enabled."""

    def test_cross_kind_fallback_to_selected(self):
        # only conda activated; container is loaded but not selected
        manager = _make_deployment_manager("conda")
        assert manager.selected_plugin_kinds == {"conda"}

        conda_spec = _make_conda_spec("env.yaml")
        container_spec = _make_container_spec("ubuntu:22.04")
        container_spec.fallback = conda_spec

        cached = MagicMock(name="cached_conda_env")
        manager._env_cache[conda_spec] = cached

        # container not selected -> falls back to conda (cached) without deploying
        assert manager.get_env(container_spec) is cached

    def test_not_selected_without_fallback_returns_none(self):
        manager = _make_deployment_manager("conda")
        container_spec = _make_container_spec("alpine:3")  # no fallback
        assert manager.get_env(container_spec) is None

    def test_selected_kind_returns_cache_directly(self):
        manager = _make_deployment_manager("conda")
        conda_spec = _make_conda_spec("env.yaml")
        cached = MagicMock(name="cached_conda_env")
        manager._env_cache[conda_spec] = cached
        # kind selected and already cached -> returns cache, no deploy attempt
        assert manager.get_env(conda_spec) is cached

    def test_multiple_kinds_can_be_activated(self):
        manager = _make_deployment_manager("conda", "container")
        assert manager.selected_plugin_kinds == {"conda", "container"}


# ---------------------------------------------------------------------------
# Integration tests — require conda executable and plugins
# ---------------------------------------------------------------------------

from .common import run, dpath, conda
from .conftest import skip_on_windows
from .common import apptainer, connected


@skip_on_windows
@conda
def test_software_directive_conda():
    """Test that the software: directive works with conda() factory."""
    run(
        dpath("test_software_directive_conda"),
        deployment_method={"conda"},
    )


@skip_on_windows
@conda
def test_software_directive_callable():
    """Test that the software: directive works with a callable envfile."""
    run(
        dpath("test_software_directive_callable"),
        deployment_method={"conda"},
    )


@skip_on_windows
@connected
@conda
def test_software_directive_within():
    """Test software: conda(envfile=..., within=container(...)) composition."""
    run(
        dpath("test_software_within"),
        deployment_method={"conda", "container"},
    )


@skip_on_windows
@conda
def test_software_directive_fallback():
    """Test software: envmodules(...) or conda(...) fallback chain.

    Both envmodules and conda plugins must be loaded so that both
    factory functions are available as globals during Snakefile parsing.
    The envmodules spec will be primary, with conda as fallback.
    On systems without Lmod, env resolution falls through to conda.
    """
    run(
        dpath("test_software_fallback"),
        deployment_method={"conda"},
    )
