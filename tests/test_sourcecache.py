from pathlib import Path
import tempfile
from unittest.mock import patch, sentinel
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
