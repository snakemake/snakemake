
Tests for handling of rule renaming and path prefixing for nested module imports.

To test locally, best create testing environment

    $ conda env create -f test-environment.yml -n snakemake-testing
    $ conda activate snakemake-testing
    $ pip install -e .

from within snakemake repository root as suggested on https://snakemake.readthedocs.io/en/stable/project_info/contributing.html#testing-guidelines and run with

    $ pip install nose
    $ nosetests tests.tests:test_module_nested


