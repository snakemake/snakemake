from pathlib import Path
import tempfile
from unittest.mock import patch
from snakemake.sourcecache import GithubFile, GitlabFile


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
