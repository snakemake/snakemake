import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

from snakemake.persistence import (
    Persistence,
    _normalize_python_code,
    _normalize_shell_code,
)


class TestCleanupContainers:
    @patch("snakemake.deployment.singularity.Image", autospec=True, create=True)
    @patch("snakemake.dag.DAG", autospec=True, create=True)
    def test_one_unrequired_container_gets_removed(self, mock_dag, mock_img):
        snakecode = """rule all:
    input:
        "foo.txt"

rule foo:
    output:
        "foo.txt"
    container:
        "docker://quay.io/mbhall88/rasusa:0.7.0"
    shell:
        "rasusa --help &> {output}"
"""
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpdirpath = Path(tmpdirname)
            singularity_dir = tmpdirpath / ".snakemake/singularity"
            singularity_dir.mkdir(parents=True)
            snakefile = tmpdirpath / "Snakefile"
            snakefile.write_text(snakecode)

            unrequired_img_path = (
                singularity_dir / "32077ccc4d05977ef8b94ee6f74073fd.simg"
            )
            unrequired_img_path.touch()
            assert unrequired_img_path.exists()

            required_img_path = (
                singularity_dir / "0581af3d3099c1cc8cf0088c8efe1439.simg"
            )
            required_img_path.touch()
            assert required_img_path.exists()
            mock_img.path = str(required_img_path)
            container_imgs = {"docker://quay.io/mbhall88/rasusa:0.7.0": mock_img}
            mock_dag.container_imgs = container_imgs
            persistence = Persistence(dag=mock_dag)
            persistence.container_img_path = str(singularity_dir)

            persistence.cleanup_containers()

            assert required_img_path.exists()
            assert not unrequired_img_path.exists()


class TestCodeChanged:
    """Tests for Persistence._code_changed with mocked dependencies."""

    def _make_persistence(self, recorded_code, fmt_version=6):
        p = MagicMock()
        p.record_format_version.return_value = fmt_version
        p.code.return_value = recorded_code
        p._code = lambda rule: Persistence._code(p, rule)
        return p

    def _make_job(self, shellcmd=None, run_func_src=None):
        job = MagicMock()
        job.rule.shellcmd = shellcmd
        job.rule.run_func_src = run_func_src
        return job

    def test_shell_rule_cosmetic_change_not_detected(self):
        recorded = "echo hello\n# old comment\n"
        current = "echo hello\n# new comment\n"
        p = self._make_persistence(recorded)
        job = self._make_job(shellcmd=current)
        assert Persistence._code_changed(p, job, file="out.txt") is False

    def test_shell_rule_real_change_detected(self):
        recorded = "echo hello\n"
        current = "echo goodbye\n"
        p = self._make_persistence(recorded)
        job = self._make_job(shellcmd=current)
        assert Persistence._code_changed(p, job, file="out.txt") is True

    def test_python_rule_cosmetic_change_not_detected(self):
        recorded = "x = 1\n# old\n"
        current = "x = 1\n# new\n"
        p = self._make_persistence(recorded)
        job = self._make_job(run_func_src=current)
        assert Persistence._code_changed(p, job, file="out.txt") is False

    def test_python_rule_real_change_detected(self):
        recorded = "x = 1\n"
        current = "x = 2\n"
        p = self._make_persistence(recorded)
        job = self._make_job(run_func_src=current)
        assert Persistence._code_changed(p, job, file="out.txt") is True

    def test_old_format_version_returns_false(self):
        p = self._make_persistence("anything", fmt_version=2)
        job = self._make_job(shellcmd="echo hi")
        assert Persistence._code_changed(p, job, file="out.txt") is False

    def test_no_format_version_returns_false(self):
        p = self._make_persistence("anything", fmt_version=None)
        job = self._make_job(shellcmd="echo hi")
        assert Persistence._code_changed(p, job, file="out.txt") is False

    def test_no_recorded_code_returns_false(self):
        p = self._make_persistence(recorded_code=None)
        job = self._make_job(shellcmd="echo hi")
        assert Persistence._code_changed(p, job, file="out.txt") is False

    def test_indented_python_rule_cosmetic_change_not_detected(self):
        """run_func_src from real Snakefiles is indented."""
        recorded = "        x = 1\n        # old\n"
        current = "        x = 1\n        # new\n"
        p = self._make_persistence(recorded)
        job = self._make_job(run_func_src=current)
        assert Persistence._code_changed(p, job, file="out.txt") is False

    def test_recorded_exists_but_current_is_none(self):
        """Script/notebook rule where _code returns None — recorded is stale."""
        p = self._make_persistence("old code")
        job = self._make_job(shellcmd=None, run_func_src=None)
        assert Persistence._code_changed(p, job, file="out.txt") is True


class TestNormalizePythonCode:
    def test_comment_change_ignored(self):
        a = "x = 1\n# old comment\ny = 2\n"
        b = "x = 1\n# new comment\ny = 2\n"
        assert _normalize_python_code(a) == _normalize_python_code(b)

    def test_whitespace_change_ignored(self):
        a = "x = 1\ny = 2\n"
        b = "x  =  1\n\n\ny  =  2\n"
        assert _normalize_python_code(a) == _normalize_python_code(b)

    def test_reformatting_ignored(self):
        a = "result = foo(a,b,c)\n"
        b = "result = foo(\n    a,\n    b,\n    c,\n)\n"
        assert _normalize_python_code(a) == _normalize_python_code(b)

    def test_variable_rename_detected(self):
        a = "x = 1\n"
        b = "y = 1\n"
        assert _normalize_python_code(a) != _normalize_python_code(b)

    def test_value_change_detected(self):
        a = "x = 1\n"
        b = "x = 2\n"
        assert _normalize_python_code(a) != _normalize_python_code(b)

    def test_logic_change_detected(self):
        a = "if x > 0:\n    print(x)\n"
        b = "if x < 0:\n    print(x)\n"
        assert _normalize_python_code(a) != _normalize_python_code(b)

    def test_indented_comment_change_ignored(self):
        """Real run_func_src is indented; normalization must handle this."""
        a = "        x = 1\n        # old comment\n        y = 2\n"
        b = "        x = 1\n        # new comment\n        y = 2\n"
        assert _normalize_python_code(a) == _normalize_python_code(b)

    def test_syntax_error_falls_back_to_raw(self):
        bad = "def foo(:\n"
        assert _normalize_python_code(bad) == bad


class TestNormalizeShellCode:
    def test_comment_lines_ignored(self):
        a = "echo hello\n# this is a comment\necho world\n"
        b = "echo hello\necho world\n"
        assert _normalize_shell_code(a) == _normalize_shell_code(b)

    def test_blank_lines_ignored(self):
        a = "echo hello\n\n\necho world\n"
        b = "echo hello\necho world\n"
        assert _normalize_shell_code(a) == _normalize_shell_code(b)

    def test_leading_trailing_whitespace_ignored(self):
        a = "  echo hello  \n  echo world  \n"
        b = "echo hello\necho world\n"
        assert _normalize_shell_code(a) == _normalize_shell_code(b)

    def test_command_change_detected(self):
        a = "echo hello\n"
        b = "echo goodbye\n"
        assert _normalize_shell_code(a) != _normalize_shell_code(b)
