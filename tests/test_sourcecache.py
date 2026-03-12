from pathlib import Path
import tempfile
from unittest.mock import patch
from snakemake.sourcecache import GithubFile, GitlabFile, infer_source_file


def test_gitlabfile_host_propagation():
    file = GitlabFile(repo="owner/repo", path="parent/path", tag="tag", host="host.com")

    assert file.get_basedir().host == file.host
    assert file.join("another/path").host == file.host


def test_github_file_fetch():
    # mock the HostedGitRepo.file_exists such that it returns False, causing GithubFile.open to trigger a fetch
    with patch("snakemake.sourcecache.HostedGitRepo.file_exists", return_value=False):
        with tempfile.TemporaryDirectory() as tmpdir:
            file = GithubFile(
                repo="snakemake/snakemake",
                path="README.md",
                branch="main",
                cache_path=Path(tmpdir),
            )
            content = file.open().read()
            assert content


def test_infer_source_file_from_github_raw_url():
    file = infer_source_file(
        "https://raw.githubusercontent.com/snakemake/snakemake/main/README.md"
    )

    assert isinstance(file, GithubFile)
    assert file.repo == "snakemake/snakemake"
    assert file.branch == "main"
    assert file.path == "README.md"


def test_infer_source_file_from_github_blob_url():
    file = infer_source_file(
        "https://github.com/snakemake/snakemake/blob/main/tests/test_multiple_includes/Snakefile"
    )

    assert isinstance(file, GithubFile)
    assert file.repo == "snakemake/snakemake"
    assert file.branch == "main"
    assert file.path == "tests/test_multiple_includes/Snakefile"


def test_infer_source_file_from_gitlab_api_url_custom_host():
    file = infer_source_file(
        "https://gitlab.example.com/api/v4/projects/group%2Fproject/repository/files/workflow%2FSnakefile/raw?ref=main"
    )

    assert isinstance(file, GitlabFile)
    assert file.host == "gitlab.example.com"
    assert file.repo == "group/project"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"
