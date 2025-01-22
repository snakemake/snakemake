import os
from pathlib import Path
import tempfile
import time
from snakemake.sourcecache import GithubFile, GitlabFile


def test_github_file():
    currtime = time.time()
    file = GithubFile(repo="snakemake/snakemake", path="README.md", tag="v8.27.1")

    with tempfile.TemporaryDirectory() as tempdir:
        file.cache_path = Path(tempdir)
        assert file.mtime() < currtime
        stored = file.get_path_or_uri()
        os.path.exists(stored)
        assert str(stored).startswith(tempdir)


def test_gitlab_file_host_propagation():
    file = GitlabFile(repo="owner/repo", path="parent/path", tag="tag", host="host.com")

    assert file.get_basedir().host == file.host
    assert file.join("another/path").host == file.host
