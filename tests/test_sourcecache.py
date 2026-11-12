import os
from pathlib import Path
import tempfile
from unittest.mock import patch, sentinel

import git
import pytest

from snakemake.exceptions import WorkflowError
from snakemake.sourcecache import (
    GenericSourceFile,
    GithubFile,
    GitlabFile,
    LocalGitFile,
    SourceCache,
    infer_source_file,
)

GIT_FILE_CONTENTS = {
    "non_utf8.bin": b"\x87\x00\xff\n",
    "final_newline.txt": b"text\n",
    "two_final_newlines.txt": b"text\n\n",
    "no_final_newline.txt": b"text",
}


@pytest.fixture
def git_repo(tmp_path):
    repo_path = tmp_path / "repo"
    repo = git.Repo.init(repo_path)
    for name, content in GIT_FILE_CONTENTS.items():
        (repo_path / name).write_bytes(content)
    repo.index.add(list(GIT_FILE_CONTENTS))
    repo.index.commit("add files")
    yield repo
    repo.close()


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


def test_infer_source_file_from_empty_host_shorthand_falls_back_to_generic():
    path = "gh::snakemake/snakemake@main"

    file = infer_source_file(path)

    assert isinstance(file, GenericSourceFile)
    assert file.get_path_or_uri(secret_free=False) == path


def test_infer_source_file_from_empty_ref_shorthand_falls_back_to_generic():
    path = "gh:snakemake/snakemake@:workflow/Snakefile"

    file = infer_source_file(path)

    assert isinstance(file, GenericSourceFile)
    assert file.get_path_or_uri(secret_free=False) == path


def test_generic_source_file_cache_path_encodes_host_port():
    file = GenericSourceFile("http://127.0.0.1:56016/Snakefile")

    assert file.get_cache_path() == os.path.join(
        "http", "127.0.0.1%3A56016", "Snakefile"
    )


def test_infer_source_file_from_github_shorthand_custom_host_without_path():
    file = infer_source_file("gh:github.example.com:snakemake/snakemake@main")

    assert isinstance(file, GithubFile)
    assert file.host == "github.example.com"
    assert file.repo == "snakemake/snakemake"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_shorthand_trailing_colon_defaults_to_snakefile():
    file = infer_source_file("gh:snakemake/snakemake@main:")

    assert isinstance(file, GithubFile)
    assert file.repo == "snakemake/snakemake"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_shorthand_no_slash_in_repo_falls_back_to_generic():
    path = "gh:notarepo@main"

    file = infer_source_file(path)

    assert isinstance(file, GenericSourceFile)
    assert file.get_path_or_uri(secret_free=False) == path


def test_github_shorthand_with_slashed_ref():
    file = infer_source_file("gh:owner/repo@perf/autobump:path/to/file.txt")

    assert isinstance(file, GithubFile)
    assert file.repo == "owner/repo"
    assert file.branch == "perf/autobump"
    assert file.path == "path/to/file.txt"


def test_gitlab_shorthand_with_slashed_ref():
    file = infer_source_file("gl:group/project@feature/branch:path/to/file.py")

    assert isinstance(file, GitlabFile)
    assert file.repo == "group/project"
    assert file.branch == "feature/branch"
    assert file.path == "path/to/file.py"


def test_hosting_provider_cache_path_keeps_slashed_ref_distinct():
    slashed_ref = GithubFile(
        repo="owner/repo",
        path="path/to/Snakefile",
        branch="perf/autobump",
    )
    plain_ref = GithubFile(
        repo="owner/repo",
        path="autobump/path/to/Snakefile",
        branch="perf",
    )

    assert slashed_ref.get_cache_path() != plain_ref.get_cache_path()
    assert "%2F" in slashed_ref.get_cache_path()


def test_hosted_repo_cache_key_distinguishes_provider_and_host():
    old_cache = GithubFile._hosted_repos
    GithubFile._hosted_repos = {}
    GitlabFile._hosted_repos = GithubFile._hosted_repos
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch(
                "snakemake.sourcecache.HostedGitRepo",
                side_effect=[
                    sentinel.github_repo,
                    sentinel.gitlab_repo,
                    sentinel.custom_repo,
                ],
            ) as hosted_repo_factory:
                github_file = GithubFile(
                    repo="owner/repo",
                    path="workflow/Snakefile",
                    branch="main",
                    host="example.com",
                    cache_path=Path(tmpdir),
                )
                gitlab_file = GitlabFile(
                    repo="owner/repo",
                    path="workflow/Snakefile",
                    branch="main",
                    host="example.com",
                    cache_path=Path(tmpdir),
                )
                custom_host_file = GithubFile(
                    repo="owner/repo",
                    path="workflow/Snakefile",
                    branch="main",
                    host="github.example.com",
                    cache_path=Path(tmpdir),
                )

                assert github_file.hosted_repo is sentinel.github_repo
                assert gitlab_file.hosted_repo is sentinel.gitlab_repo
                assert custom_host_file.hosted_repo is sentinel.custom_repo
                assert hosted_repo_factory.call_count == 3
    finally:
        GithubFile._hosted_repos = old_cache
        GitlabFile._hosted_repos = old_cache


@pytest.mark.parametrize("name", GIT_FILE_CONTENTS)
def test_local_git_file_preserves_bytes(git_repo, tmp_path, name):
    file = LocalGitFile(git_repo.working_dir, name, commit=git_repo.head.commit.hexsha)

    with SourceCache(tmp_path / "cache").open(file, "rb") as content:
        assert content.read() == GIT_FILE_CONTENTS[name]


def make_hosted_file(git_repo, cache_path, path):
    """Return a GithubFile served from a local bare clone of git_repo.

    HostedGitRepo reuses an existing clone instead of cloning from the
    network, so placing one where it looks makes the test run offline.
    """
    file = GithubFile(repo="owner/repo", path=path, commit=git_repo.head.commit.hexsha)
    file.cache_path = cache_path
    git.Repo.clone_from(
        git_repo.working_dir, file.cache_path / file.host / file.repo, bare=True
    ).close()
    return file


@pytest.mark.parametrize("name", GIT_FILE_CONTENTS)
def test_hosted_git_file_preserves_bytes(git_repo, tmp_path, name):
    source_cache = SourceCache(tmp_path / "cache")
    file = make_hosted_file(git_repo, source_cache.cache_path, name)

    with (
        patch.dict(GithubFile._hosted_repos, clear=True),
        source_cache.open(file, "rb") as content,
    ):
        assert content.read() == GIT_FILE_CONTENTS[name]


def test_local_git_file_missing_path_raises_workflow_error(git_repo, tmp_path):
    file = LocalGitFile(
        git_repo.working_dir, "missing.txt", commit=git_repo.head.commit.hexsha
    )

    with pytest.raises(WorkflowError, match="Failed to get local git source file"):
        SourceCache(tmp_path / "cache").open(file, "rb")


def test_hosted_git_file_missing_path_raises_workflow_error(git_repo, tmp_path):
    file = make_hosted_file(git_repo, tmp_path / "cache", "missing.txt")

    with (
        patch.dict(GithubFile._hosted_repos, clear=True),
        pytest.raises(WorkflowError, match="Failed to get git source file"),
    ):
        file.open()
