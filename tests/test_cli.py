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


# ---------------------------------------------------------------------------
# Legacy fallthrough
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Subcommand routing
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Lint subcommand
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Info subcommand group
# ---------------------------------------------------------------------------


def test_info_rules(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.info.run_list_rules") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["info", "rules"])

        assert result.exit_code == 0
        mock_op.assert_called_once()
        # should not pass only_targets
        assert mock_op.call_args.kwargs.get("only_targets", False) is False


def test_info_target_rules(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.info.run_list_rules") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["info", "target-rules"])

        assert result.exit_code == 0
        mock_op.assert_called_once()
        assert mock_op.call_args.kwargs.get("only_targets", False) is True


def test_info_summary(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.info.run_summary") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["info", "summary"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag, detailed=False)


def test_info_summary_detailed(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.info.run_summary") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["info", "summary", "--detailed"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag, detailed=True)


def test_info_changes_code(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.info.run_list_changes") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        from snakemake.settings.enums import ChangeType

        result = runner.invoke(cli, ["info", "changes", "code"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag, ChangeType.CODE)


def test_info_changes_input(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.info.run_list_changes") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        from snakemake.settings.enums import ChangeType

        result = runner.invoke(cli, ["info", "changes", "input"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag, ChangeType.INPUT)


def test_info_changes_params(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.info.run_list_changes") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        from snakemake.settings.enums import ChangeType

        result = runner.invoke(cli, ["info", "changes", "params"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag, ChangeType.PARAMS)


def test_info_changes_invalid_type(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with patch("snakemake.cli.common.SnakemakeApi", mock_cls):
        result = runner.invoke(cli, ["info", "changes", "invalid"])

        assert result.exit_code != 0


def test_info_untracked(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.info.run_list_untracked") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["info", "untracked"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag)


def test_info_compilation(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.info.run_print_compilation") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["info", "compilation"])

        assert result.exit_code == 0
        mock_op.assert_called_once()


def test_info_with_snakefile(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.info.run_list_rules"),
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["info", "rules", "-s", "workflow/Snakefile"])

        assert result.exit_code == 0
        assert mock_api.workflow.call_args.kwargs["snakefile"] == Path(
            "workflow/Snakefile"
        )


# ---------------------------------------------------------------------------
# DAG subcommand
# ---------------------------------------------------------------------------


def test_dag_default(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.dagviz.run_printdag") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["dagviz"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag)


def test_dag_rulegraph(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.dagviz.run_printrulegraph") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["dagviz", "--rulegraph"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag)


def test_dag_filegraph(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.dagviz.run_printfilegraph") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["dagviz", "--filegraph"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag)


def test_dag_d3(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.dagviz.run_printd3dag") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["dagviz", "--d3"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag)


def test_dag_with_snakefile(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.dagviz.run_printdag"),
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["dagviz", "-s", "workflow/Snakefile"])

        assert result.exit_code == 0
        assert mock_api.workflow.call_args.kwargs["snakefile"] == Path(
            "workflow/Snakefile"
        )


# ---------------------------------------------------------------------------
# Clean subcommand
# ---------------------------------------------------------------------------


def test_clean_all(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.clean.run_delete_output") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["clean", "--all"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag, dryrun=False)


def test_clean_temp(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.clean.run_delete_output") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["clean", "--temp"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag, only_temp=True, dryrun=False)


def test_clean_dry_run(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.clean.run_delete_output") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["clean", "--all", "--dry-run"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag, dryrun=True)


def test_clean_shadow(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.clean.run_cleanup_shadow") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["clean", "--shadow"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag)


def test_clean_metadata(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.clean.run_cleanup_metadata") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(
            cli, ["clean", "--metadata", "file1.txt", "--metadata", "file2.txt"]
        )

        assert result.exit_code == 0
        mock_op.assert_called_once_with(
            mock_dag, [Path("file1.txt"), Path("file2.txt")]
        )


def test_clean_no_flags_errors(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with patch("snakemake.cli.common.SnakemakeApi", mock_cls):
        result = runner.invoke(cli, ["clean"])

        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# Utils subcommand group
# ---------------------------------------------------------------------------


def test_utils_containerize(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.utils.run_containerize") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["utils", "containerize"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag)


def test_utils_generate_unit_tests_default_path(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.utils.run_generate_unit_tests") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["utils", "generate-unit-tests"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag, Path(".tests/unit"))


def test_utils_generate_unit_tests_custom_path(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.utils.run_generate_unit_tests") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["utils", "generate-unit-tests", "my/tests"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag, Path("my/tests"))


def test_utils_archive(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    mock_workflow = MagicMock()
    mock_dag = MagicMock()
    mock_api.workflow.return_value = mock_workflow
    mock_workflow.dag.return_value = mock_dag

    with (
        patch("snakemake.cli.commands.utils.run_archive") as mock_op,
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(cli, ["utils", "archive", "workflow.tar.gz"])

        assert result.exit_code == 0
        mock_op.assert_called_once_with(mock_dag, Path("workflow.tar.gz"))


def test_utils_archive_requires_file_arg(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with patch("snakemake.cli.common.SnakemakeApi", mock_cls):
        result = runner.invoke(cli, ["utils", "archive"])

        assert result.exit_code != 0


def test_utils_containerize_with_snakefile(runner):
    mock_cls, mock_api = _mock_snakemake_api()
    with (
        patch("snakemake.cli.commands.utils.run_containerize"),
        patch("snakemake.cli.common.SnakemakeApi", mock_cls),
    ):
        result = runner.invoke(
            cli, ["utils", "containerize", "-s", "workflow/Snakefile"]
        )

        assert result.exit_code == 0
        assert mock_api.workflow.call_args.kwargs["snakefile"] == Path(
            "workflow/Snakefile"
        )


# ---------------------------------------------------------------------------
# Help text
# ---------------------------------------------------------------------------


def test_lint_help(runner):
    result = runner.invoke(cli, ["lint", "--help"])
    assert result.exit_code == 0
    assert "--snakefile" in result.output
    assert "--directory" in result.output


def test_unlock_help(runner):
    result = runner.invoke(cli, ["unlock", "--help"])
    assert result.exit_code == 0
    assert "--snakefile" in result.output
    assert "--directory" in result.output


def test_info_help(runner):
    result = runner.invoke(cli, ["info", "--help"])
    assert result.exit_code == 0
    assert "rules" in result.output
    assert "summary" in result.output
    assert "changes" in result.output


def test_dag_help(runner):
    result = runner.invoke(cli, ["dagviz", "--help"])
    assert result.exit_code == 0
    assert "--rulegraph" in result.output
    assert "--filegraph" in result.output


def test_clean_help(runner):
    result = runner.invoke(cli, ["clean", "--help"])
    assert result.exit_code == 0
    assert "--all" in result.output
    assert "--temp" in result.output
    assert "--shadow" in result.output
    assert "--metadata" in result.output


def test_utils_help(runner):
    result = runner.invoke(cli, ["utils", "--help"])
    assert result.exit_code == 0
    assert "containerize" in result.output
    assert "generate-unit-tests" in result.output
    assert "archive" in result.output