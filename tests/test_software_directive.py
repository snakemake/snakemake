__authors__ = ["Yun Jiang"]
__copyright__ = "Copyright 2026, Yun Jiang"
__license__ = "MIT"

import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(__file__))

from snakemake_interface_software_deployment_plugins import EnvSpecSourceFile
from snakemake_software_deployment_plugin_conda import EnvSpec as CondaEnvSpec
from snakemake_software_deployment_plugin_container import EnvSpec as ContainerEnvSpec
from snakemake_software_deployment_plugin_envmodules import EnvSpec as EnvModuleEnvSpec

from snakemake.deployment import EnvSpecs


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
@apptainer
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
        deployment_method={"envmodules", "conda"},
    )
