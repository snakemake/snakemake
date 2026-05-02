from pathlib import Path
import tempfile
from unittest.mock import patch, MagicMock
from snakemake.sourcecache import (
    GenericSourceFile,
    GithubFile,
    GitlabFile,
    _build_ref_path_candidates,
    _infer_github_file,
    _infer_gitlab_file,
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


def test_infer_source_file_from_github_shorthand_custom_host_without_path():
    file = infer_source_file("gh:github.example.com:snakemake/snakemake@main")

    assert isinstance(file, GithubFile)
    assert file.host == "github.example.com"
    assert file.repo == "snakemake/snakemake"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_shorthand_trailing_colon_defaults_to_snakefile():
    # gh:repo@ref: (trailing colon, empty path) should default to workflow/Snakefile
    file = infer_source_file("gh:snakemake/snakemake@main:")

    assert isinstance(file, GithubFile)
    assert file.repo == "snakemake/snakemake"
    assert file.branch == "main"
    assert file.path == "workflow/Snakefile"


def test_infer_source_file_from_shorthand_no_slash_in_repo_falls_back_to_generic():
    # gh:notarepo@main has no slash in the repo part — should fall back gracefully
    path = "gh:notarepo@main"

    file = infer_source_file(path)

    assert isinstance(file, GenericSourceFile)
    assert file.get_path_or_uri(secret_free=False) == path


# ---------------------------------------------------------------------------
# Tests for ambiguous branch/path parsing (branch names containing slashes)
# ---------------------------------------------------------------------------


def test_build_ref_path_candidates():
    # Simulates path_parts for: /owner/repo/raw/perf/autobump/.editorconfig
    path_parts = ["owner", "repo", "raw", "perf", "autobump", ".editorconfig"]
    candidates = _build_ref_path_candidates(path_parts, ref_start=3)

    # 2 candidates: the loop leaves at least one segment for the path
    assert len(candidates) == 2
    assert candidates[0] == ("perf", "autobump/.editorconfig")
    assert candidates[1] == ("perf/autobump", ".editorconfig")


def test_build_ref_path_candidates_single_segment():
    # Only one segment after ref_start → only one candidate, no ambiguity.
    path_parts = ["owner", "repo", "raw", "main", "file.txt"]
    candidates = _build_ref_path_candidates(path_parts, ref_start=3)

    assert len(candidates) == 1
    assert candidates[0] == ("main", "file.txt")


def test_infer_github_file_ambiguous_branch_creates_candidates():
    """GitHub URL with a branch containing a slash should store candidates."""
    file = _infer_github_file(
        "https://github.com/snakemake-workflows/dna-seq-varlociraptor/raw/perf/autobump/.editorconfig"
    )

    assert isinstance(file, GithubFile)
    assert file.repo == "snakemake-workflows/dna-seq-varlociraptor"
    # Naive split takes the first segment as the branch
    assert file.branch == "perf"
    assert file.path == "autobump/.editorconfig"
    # But candidates are stored for lazy resolution
    assert file._ref_path_candidates is not None
    assert len(file._ref_path_candidates) == 2
    assert ("perf", "autobump/.editorconfig") in file._ref_path_candidates
    assert ("perf/autobump", ".editorconfig") in file._ref_path_candidates


def test_infer_github_file_raw_githubusercontent_ambiguous_branch():
    """raw.githubusercontent.com URL with a slashed branch creates candidates."""
    file = _infer_github_file(
        "https://raw.githubusercontent.com/owner/repo/perf/autobump/.editorconfig"
    )

    assert isinstance(file, GithubFile)
    assert file.repo == "owner/repo"
    assert file.branch == "perf"
    assert file.path == "autobump/.editorconfig"
    assert file._ref_path_candidates is not None
    assert ("perf/autobump", ".editorconfig") in file._ref_path_candidates


def test_infer_github_file_no_ambiguity_no_candidates():
    """A simple branch name should not create candidates."""
    file = _infer_github_file("https://github.com/owner/repo/raw/main/README.md")

    assert isinstance(file, GithubFile)
    assert file.branch == "main"
    assert file.path == "README.md"
    assert file._ref_path_candidates is None


def test_infer_github_file_encoded_slash_no_candidates():
    """URL-encoded slash in branch → single path segment, no ambiguity."""
    file = _infer_github_file(
        "https://github.com/owner/repo/blob/feature%2Ffoo/workflow/Snakefile"
    )

    assert isinstance(file, GithubFile)
    assert file.branch == "feature/foo"
    assert file.path == "workflow/Snakefile"
    # URL-encoded slash is kept together → only one candidate per segment
    # 3 segments after "blob": "feature%2Ffoo", "workflow", "Snakefile"
    # candidates = [("feature/foo", "workflow/Snakefile"), ("feature/foo/workflow", "Snakefile")]
    # → 2 candidates, but the first one is already correct anyway
    assert file._ref_path_candidates is not None
    assert file._ref_path_candidates[0] == ("feature/foo", "workflow/Snakefile")


def test_infer_gitlab_file_ambiguous_branch_creates_candidates():
    """GitLab URL with a slashed branch should store candidates."""
    file = _infer_gitlab_file(
        "https://gitlab.example.com/group/project/-/raw/perf/autobump/.editorconfig"
    )

    assert isinstance(file, GitlabFile)
    assert file.repo == "group/project"
    assert file.host == "gitlab.example.com"
    assert file.branch == "perf"
    assert file.path == "autobump/.editorconfig"
    assert file._ref_path_candidates is not None
    assert ("perf/autobump", ".editorconfig") in file._ref_path_candidates


def test_infer_gitlab_file_no_ambiguity_no_candidates():
    """A URL with only one path segment after the branch is unambiguous."""
    file = _infer_gitlab_file(
        "https://gitlab.example.com/group/project/-/blob/main/Snakefile"
    )

    assert isinstance(file, GitlabFile)
    assert file.branch == "main"
    assert file.path == "Snakefile"
    assert file._ref_path_candidates is None


def test_resolve_ref_picks_correct_candidate():
    """_resolve_ref should update branch and path to the matching candidate."""
    file = _infer_github_file(
        "https://github.com/owner/repo/raw/perf/autobump/.editorconfig"
    )

    assert file._ref_path_candidates is not None

    # Mock hosted_repo.ref_exists: "perf" does not exist, "perf/autobump" does
    mock_hosted_repo = MagicMock()
    mock_hosted_repo.ref_exists.side_effect = lambda commit, branch_or_tag: (
        branch_or_tag == "perf/autobump"
    )

    with patch.object(
        type(file),
        "hosted_repo",
        new_callable=lambda: property(lambda self: mock_hosted_repo),
    ):
        file._resolve_ref()

    assert file.branch == "perf/autobump"
    assert file.path == ".editorconfig"
    assert file._ref_path_candidates is None


def test_resolve_ref_keeps_naive_split_when_no_candidate_matches():
    """If no candidate ref exists, _resolve_ref keeps the naive (first-segment) split."""
    file = _infer_github_file(
        "https://github.com/owner/repo/raw/perf/autobump/.editorconfig"
    )

    mock_hosted_repo = MagicMock()
    mock_hosted_repo.ref_exists.return_value = False

    with patch.object(
        type(file),
        "hosted_repo",
        new_callable=lambda: property(lambda self: mock_hosted_repo),
    ):
        file._resolve_ref()

    # Keeps the naive split
    assert file.branch == "perf"
    assert file.path == "autobump/.editorconfig"
    assert file._ref_path_candidates is None


def test_resolve_ref_gitlab_picks_correct_candidate():
    """_resolve_ref works for GitLab files too."""
    file = _infer_gitlab_file(
        "https://gitlab.example.com/group/project/-/raw/feature/branch/path/to/file.py"
    )

    assert file._ref_path_candidates is not None

    mock_hosted_repo = MagicMock()
    mock_hosted_repo.ref_exists.side_effect = lambda commit, branch_or_tag: (
        branch_or_tag == "feature/branch"
    )

    with patch.object(
        type(file),
        "hosted_repo",
        new_callable=lambda: property(lambda self: mock_hosted_repo),
    ):
        file._resolve_ref()

    assert file.branch == "feature/branch"
    assert file.path == "path/to/file.py"


# ---------------------------------------------------------------------------
# Tests for get_path_or_uri round-tripping
# ---------------------------------------------------------------------------


def test_github_file_get_path_or_uri_round_trip():
    """GithubFile.get_path_or_uri should produce a valid URL containing repo, ref, and path."""
    file = GithubFile(repo="owner/repo", path="workflow/Snakefile", branch="main")
    uri = file.get_path_or_uri(secret_free=True)
    assert "owner/repo" in uri
    assert "main" in uri
    assert "workflow/Snakefile" in uri


def test_gitlab_file_get_path_or_uri_round_trip():
    """GitlabFile.get_path_or_uri should produce a valid API URL."""
    file = GitlabFile(
        repo="group/project",
        path="workflow/Snakefile",
        branch="main",
        host="gitlab.example.com",
    )
    uri = file.get_path_or_uri(secret_free=True)
    assert "gitlab.example.com" in uri
    assert "ref=main" in uri


def test_github_shorthand_with_slashed_ref():
    """Shorthand notation unambiguously handles refs with slashes."""
    file = infer_source_file("gh:owner/repo@perf/autobump:path/to/file.txt")

    assert isinstance(file, GithubFile)
    assert file.repo == "owner/repo"
    assert file.branch == "perf/autobump"
    assert file.path == "path/to/file.txt"
    # Shorthand uses : as delimiter — no ambiguity, no candidates
    assert file._ref_path_candidates is None


def test_gitlab_shorthand_with_slashed_ref():
    """Shorthand notation unambiguously handles refs with slashes for GitLab."""
    file = infer_source_file("gl:group/project@feature/branch:path/to/file.py")

    assert isinstance(file, GitlabFile)
    assert file.repo == "group/project"
    assert file.branch == "feature/branch"
    assert file.path == "path/to/file.py"
    assert file._ref_path_candidates is None
