from pathlib import Path
import tempfile
from unittest.mock import patch
from snakemake.sourcecache import (
    GenericSourceFile,
    GithubFile,
    GitlabFile,
    infer_source_file,
)


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


def test_infer_source_file_from_github_raw_url_on_github_com():
    file = infer_source_file(
        "https://github.com/snakemake/snakemake/raw/main/README.md"
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


def test_infer_source_file_from_github_blob_url_with_encoded_slash_ref():
    file = infer_source_file(
        "https://github.com/snakemake/snakemake/blob/feature%2Ffoo/workflow/Snakefile"
    )

    assert isinstance(file, GithubFile)
    assert file.repo == "snakemake/snakemake"
    assert file.branch == "feature/foo"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_gitlab_api_url_custom_host():
    file = infer_source_file(
        "https://gitlab.example.com/api/v4/projects/group%2Fproject/repository/files/workflow%2FSnakefile/raw?ref=main"
    )

    assert isinstance(file, GitlabFile)
    assert file.host == "gitlab.example.com"
    assert file.repo == "group/project"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_gitlab_blob_url_custom_host():
    file = infer_source_file(
        "https://gitlab.example.com/group/project/-/blob/main/workflow/Snakefile"
    )

    assert isinstance(file, GitlabFile)
    assert file.host == "gitlab.example.com"
    assert file.repo == "group/project"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_gitlab_raw_url_custom_host():
    file = infer_source_file(
        "https://gitlab.example.com/group/subgroup/project/-/raw/main/workflow/Snakefile"
    )

    assert isinstance(file, GitlabFile)
    assert file.host == "gitlab.example.com"
    assert file.repo == "group/subgroup/project"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_gitlab_blob_url_with_encoded_slash_ref():
    file = infer_source_file(
        "https://gitlab.example.com/group/project/-/blob/feature%2Ffoo/workflow/Snakefile"
    )

    assert isinstance(file, GitlabFile)
    assert file.host == "gitlab.example.com"
    assert file.repo == "group/project"
    assert file.branch == "feature/foo"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_gitlab_api_url_with_private_token_falls_back_to_generic():
    url = (
        "https://gitlab.example.com/api/v4/projects/group%2Fproject/repository/files/"
        "workflow%2FSnakefile/raw?ref=main&private_token=secret"
    )
    file = infer_source_file(url)

    assert isinstance(file, GenericSourceFile)
    assert file.get_path_or_uri(secret_free=False) == url


def test_infer_source_file_from_github_shorthand():
    file = infer_source_file("gh:snakemake/snakemake@main:workflow/Snakefile")

    assert isinstance(file, GithubFile)
    assert file.host == "github.com"
    assert file.repo == "snakemake/snakemake"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_github_shorthand_without_path():
    file = infer_source_file("gh:snakemake/snakemake@main")

    assert isinstance(file, GithubFile)
    assert file.host == "github.com"
    assert file.repo == "snakemake/snakemake"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_github_shorthand_custom_host():
    file = infer_source_file(
        "gh:github.example.com:snakemake/snakemake@main:workflow/Snakefile"
    )

    assert isinstance(file, GithubFile)
    assert file.host == "github.example.com"
    assert file.repo == "snakemake/snakemake"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"
    assert "github.example.com" in file.get_path_or_uri(secret_free=True)


def test_infer_source_file_from_gitlab_shorthand():
    file = infer_source_file("gl:group/project@main:workflow/Snakefile")

    assert isinstance(file, GitlabFile)
    assert file.host == "gitlab.com"
    assert file.repo == "group/project"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_gitlab_shorthand_without_path():
    file = infer_source_file("gl:group/project@main")

    assert isinstance(file, GitlabFile)
    assert file.host == "gitlab.com"
    assert file.repo == "group/project"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_gitlab_shorthand_custom_host():
    file = infer_source_file(
        "gl:gitlab.cern.ch:group/subgroup/project@main:workflow/Snakefile"
    )

    assert isinstance(file, GitlabFile)
    assert file.host == "gitlab.cern.ch"
    assert file.repo == "group/subgroup/project"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_old_host_first_shorthand_falls_back_to_generic():
    path = "github.com:snakemake/snakemake@main:workflow/Snakefile"

    file = infer_source_file(path)

    assert isinstance(file, GenericSourceFile)
    assert file.get_path_or_uri(secret_free=False) == path
