# Unit and regression tests for Snakemake

## Running the tests

All tests are written to by run by [pytest](https://docs.pytest.org/en/stable/). If you simply
run `pytest` it will scan for and attempt to run every test in the directory, but you will almost
certainly get errors as not all tests are working and current. See the section below.

Assuming you are on a Linux system, and have a working Conda installation, you can set up to
run the tests like so:

1. Ensure your Conda config (~/.condarc) has these options:

```
channels:
  - conda-forge
  - bioconda
channel_priority: strict
solver: libmamba
```

*note: if you need to keep different settings in your personal `~/.condarc` for some reason,
you can make the new env, activate it, and put the above into `$CONDA_PREFIX/.condarc`
so it will apply to just that environment.*

2. After checking out the branch you want to test, run these commands:

```
$ conda env create -n snakemake_test_env -f test-environment.yml
$ conda activate snakemake_test_env
$ pip install -e .
```

You may want to set a specific Python version by editing the constraint in `test-environment.yml`
before doing this.

3. You can now try running a single specific test, eg:

```
$ pytest tests/tests.py::test_modules_all
```

You can also use the `-k` flag to select tests by substring match, rather than by the full name,
and the `--co` option to preview which tests will be run. Try, for example:

```
$ pytest --co tests/tests.py -k test_modules_all
```

## Running the full test suite

The core test suite is the set of tests run as a GitHub action by the code under
`.github/workflows/main.yml`, so you should look in this file for the list of tests actually
expected to pass in a regular test environment. At the time of writing this text, the suite is:

```
tests/tests.py
tests/tests_using_conda.py
tests/test_expand.py
tests/test_io.py
tests/test_schema.py
tests/test_linting.py
tests/test_executor_test_suite.py
tests/test_api.py
tests/test_internals.py
```

Other tests in the directory may or may not work.

## Warnings and oddities

You will likely see warnings related to deprecated functions in dependent libraries, especially
botocore.

You may also get intermittent failures from tests that rely on external connectivity. The default
test suite makes connections to multiple external services.

Tests that require singularity will be auto-skipped if no singularity or apptainer installation is
available. At the time of writing neither the `singularity` package on conda-forge nor the
`apptainer` package are reliable, in that there are multiple failing tests on a standard Ubuntu
system. This is likely due to system security profiles that conda, being a non-root application,
cannot change.
The Debian/Ubuntu `singularity-container` DEB package, which must be installed by the system
administrator, does work. The equivalent RPM package should also work on RedHat-type systems.

Tests in `tests/test_api.py` require a working `git`. This is not included in
`test-environment.yml` as it's assumed you must have GIT installed to be working on the source
code, but installing git into the conda environment should work if need be.

Depending on how the Snakemake code was downloaded and installed in the test environmant,
Snakemake may not be able to determine its own version and may think that it is version 0.
The unit tests should all cope with this but be aware that you can not rely on anything that
depends on version checks.

