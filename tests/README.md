# Unit and regression tests for Snakemake

## Running the tests

All tests are written to by run by [pytest](https://docs.pytest.org/en/stable/). If you simply
run `pytest` it will discover and attempt to run every test in the directory, but you will almost
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

2. After checkoing out the branch you want to tes, run these commands:

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

The core test suite is run as a GitHub action by the code under `.github/workflows/main.yml`, so
you should look in this file for the list of tests actually expected to pass in a regular test
environment. At the time of writing this suite is:

```
past-y
```

If you want to run the entire suite, you will need to:

TODO

## Warnings and oddities

-----

Earlier notes...

I wanted to diagnose an integration test failure. Specifically, one in
tests/tests.py::test_modules_all picked up by a GitHub test runner.

I'd like to run this test in isolation on my own system (Ubuntu Linux with a working Conda-Forge
setup) in order to be able to attempt a fix for the error.

I make a suitable conda env:

```
$ cd ~/git_workspace/snakemake
$ conda env create -n snakemake8_test -f test-environment.yml
$ conda activate snakemake8_test
$ pip install -e .
```

Now how do I run this specific test? What if I ask to discover all tests?

```
$ pytest --co
```

I get a couple of errors. I think I can ignore these and specify a single test file.

```
$ pytest --co tests/tests.py -k test_modules_all
```

This is happy but it picks up three tests that have this within the name, not just the one test.

```
$ pytest --co tests/tests.py::test_modules_all
```

I struggled to find the syntax above, but this is the one. I can now run the single test. :-)

And I get the same error as in GIT. Progress!

Having fixed this, I try running all the tests. I get a bunch of failures and a bunch of warnings:

================================= short test summary info =====================================================
FAILED tests/tests.py::test_set_resources_complex_profile - AssertionError: expected successful execution
FAILED tests/tests.py::test_benchmark - AssertionError: expected successful execution
FAILED tests/tests.py::test_benchmark_jsonl - AssertionError: expected successful execution
FAILED tests/tests.py::test_script_rs - AssertionError: expected successful execution
FAILED tests/tests.py::test_profile - AssertionError: expected successful execution
FAILED tests/tests.py::test_module_complex - AssertionError: expected successful execution
FAILED tests/tests.py::test_module_complex2 - AssertionError: expected successful execution
FAILED tests/tests.py::test_module_with_script - AssertionError: expected successful execution
FAILED tests/tests.py::test_strict_mode - AssertionError: expected error on execution
FAILED tests/tests.py::test_workflow_profile - AssertionError: expected successful execution
FAILED tests/tests.py::test_no_workflow_profile - AssertionError: expected successful execution
FAILED tests/tests.py::test_runtime_conversion_from_workflow_profile - AssertionError: expected successful execution
FAILED tests/tests.py::test_resource_string_in_cli_or_profile - AssertionError: expected successful execution
FAILED tests/tests.py::test_shell_exec - AssertionError: expected successful execution
=========== 14 failed, 299 passed, 6 skipped, 525 warnings in 1419.75s (0:23:39) ==============================

I think some or all of these are because I do not have a working Singularity, but I think a major
cleanup of the tests would be a good goal for the Hackathon.

Likely:

1) The info above should be in the README.md for the tests/ subdir.
2) The tests that need singularity should be split out, and should fail fast if singularity
   is not working at all.
3) The 'import PIL' issue should be addressed, and the other one, so that
   'pytest --co' should work, even if the tests then fail due to import errors.
4) We should have a boilerplate for adding unit tests, like what I have for my own code
5) The various deprecation warnings should be addressed now before they bite us
6) Should run as many tests as possible using only conda-forge packages, in particular:
  - Replace 'stress' with 'stress-ng'
  - Avoid use of wget or wrap it somehow
  - Also consider the need for git, openmpi-bin, libopenmpi-dev, apptainer
7) Avoid any tests changing files that are comitted in GIT. I need to find out which specific
   test is changing tests/test_prebuilt_conda_script/dummy_package/src/dummy.egg-info/PKG-INFO
   and sort that one. Also, what is writing ./version.txt


## Tidying up the tests

There is a lot to do, so this might be quite scattergun. A massive amount of Yak shaving is
going to be needed. I think what I'll do is work in a new branch based off of remove-datrie as
that is where I'm starting. I'll fix what I can, then go back and re-make all the changes in a
logical-looking order, to make a PR. And then clean up this README.
