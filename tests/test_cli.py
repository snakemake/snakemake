from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from snakemake.cli.click import cli


@pytest.fixture
def runner():
    return CliRunner()


def _mock_snakemake_api():
    mock_api = MagicMock()
    mock_cls = MagicMock()
    mock_cls.return_value.__enter__ = MagicMock(return_value=mock_api)
    mock_cls.return_value.__exit__ = MagicMock(return_value=False)
    return mock_cls, mock_api


def test_no_args_falls_through(runner):
    with patch("snakemake.cli.click.legacy_main") as mock_legacy:
        mock_legacy.return_value = None
        runner.invoke(cli, [])
        mock_legacy.assert_called_once()


def test_unknown_flags_fall_through(runner):
    with patch("snakemake.cli.click.legacy_main") as mock_legacy:
        mock_legacy.return_value = None
        runner.invoke(cli, ["--dry-run"])
        mock_legacy.assert_called_once()


def test_positional_targets_fall_through(runner):
    with patch("snakemake.cli.click.legacy_main") as mock_legacy:
        mock_legacy.return_value = None
        runner.invoke(cli, ["all"])
        mock_legacy.assert_called_once()


def test_legacy_lint_flag_falls_through(runner):
    with patch("snakemake.cli.click.legacy_main") as mock_legacy:
        mock_legacy.return_value = None
        runner.invoke(cli, ["--lint"])
        mock_legacy.assert_called_once()


def test_legacy_unlock_flag_falls_through(runner):
    with patch("snakemake.cli.click.legacy_main") as mock_legacy:
        mock_legacy.return_value = None
        runner.invoke(cli, ["--unlock"])
        mock_legacy.assert_called_once()


def test_complex_legacy_invocation_falls_through(runner):
    with patch("snakemake.cli.click.legacy_main") as mock_legacy:
        mock_legacy.return_value = None
        runner.invoke(
            cli, ["--cores", "4", "--dry-run", "--snakefile", "workflow/Snakefile"]
        )
        mock_legacy.assert_called_once()


def test_lint_is_subcommand_not_target(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.click.legacy_main") as mock_legacy,
        patch("snakemake.cli.commands.lint.run_lint", return_value=True),
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        runner.invoke(cli, ["lint"])
        mock_legacy.assert_not_called()


def test_unlock_is_subcommand_not_target(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.click.legacy_main") as mock_legacy,
        patch("snakemake.cli.commands.unlock.run_unlock"),
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        runner.invoke(cli, ["unlock"])
        mock_legacy.assert_not_called()


def test_lint_default(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.lint.run_lint", return_value=True),
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["lint"])

        assert result.exit_code == 0
        call_kwargs = mock_api.workflow.call_args.kwargs
        assert call_kwargs["snakefile"] is None
        assert call_kwargs["workdir"] is None


def test_lint_with_snakefile(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.lint.run_lint", return_value=True),
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["lint", "-s", "workflow/Snakefile"])

        assert result.exit_code == 0
        assert mock_api.workflow.call_args.kwargs["snakefile"] == Path(
            "workflow/Snakefile"
        )


def test_lint_with_directory(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.lint.run_lint", return_value=True),
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["lint", "-d", "/some/workdir"])

        assert result.exit_code == 0
        assert mock_api.workflow.call_args.kwargs["workdir"] == Path("/some/workdir")


def test_lint_failure_exits_nonzero(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.lint.run_lint", return_value=False),
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["lint"])

        assert result.exit_code == 1


# ---------------------------------------------------------------------------
# Unlock subcommand
# ---------------------------------------------------------------------------


def test_unlock_default(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.unlock.run_unlock") as mock_run_unlock,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["unlock"])

        assert result.exit_code == 0
        mock_api.workflow.assert_called_once()
        mock_run_unlock.assert_called_once()


def test_unlock_with_directory(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.unlock.run_unlock"),
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["unlock", "-d", "/some/workdir"])

        assert result.exit_code == 0
        assert mock_api.workflow.call_args.kwargs["workdir"] == Path("/some/workdir")


def test_unlock_with_snakefile(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.unlock.run_unlock"),
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["unlock", "-s", "workflow/Snakefile"])

        assert result.exit_code == 0
        assert mock_api.workflow.call_args.kwargs["snakefile"] == Path(
            "workflow/Snakefile"
        )


def test_unlock_calls_dag_api(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow_api = MagicMock()
    mock_dag_api = MagicMock()
    mock_api.workflow.return_value = mock_workflow_api
    mock_workflow_api.dag.return_value = mock_dag_api

    with (
        patch("snakemake.cli.commands.unlock.run_unlock") as mock_run_unlock,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        runner.invoke(cli, ["unlock"])

        mock_workflow_api.dag.assert_called_once()
        mock_run_unlock.assert_called_once_with(mock_dag_api)
