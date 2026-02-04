"""
conftest.py for unit testing of rules generated with Snakemake 9.8.2.dev50.
"""

from pytest import fixture


def pytest_addoption(parser):
    parser.addoption("--conda-prefix", action="store", default=None)


@fixture()
def conda_prefix(request):
    conda_prefix = request.config.option.conda_prefix
    if conda_prefix:
        return ["--conda-prefix", conda_prefix]
    else:
        return []
