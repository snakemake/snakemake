from snakemake.sourcecache import GitlabFile, HttpFile


def test_GitlabFile_host_propagation():
    file = GitlabFile(repo="owner/repo", path="parent/path", tag="tag", host="host.com")

    assert file.get_basedir().host == file.host
    assert file.join("another/path").host == file.host


def test_http_file():
    file = HttpFile("http://ipv4.download.thinkbroadband.com:8080/5MB.zip")

    assert file.mtime()
