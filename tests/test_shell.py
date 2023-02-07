from snakemake.shell import shell


def test_collect_stdout():
    output = shell("echo test output", read=True)
    assert output == b"test output\n"
