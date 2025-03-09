.. _project_info-contributing:

============
Contributing
============

Contributions are welcome, and they are greatly appreciated!
Every little bit helps, and credit will always be given.

You can contribute in many ways:


----------------------
Types of Contributions
----------------------


Report Bugs
===========

Report bugs at https://github.com/snakemake/snakemake/issues

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.


Fix Bugs
========

Look through the Github issues for bugs.
If you want to start working on a bug then please write short message on the issue tracker to prevent duplicate work.


Implement Features
==================

Look through the Github issues for features.
If you want to start working on an issue then please write short message on the issue tracker to prevent duplicate work.

Contributing a plugin
~~~~~~~~~~~~~~~~~~~~~

Currently, Snakemake supports executor plugins and storage plugins.
The `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog>`_ shows which plugins are available and how to contribute new ones.

Write Documentation
===================

Snakemake could always use more documentation, whether as part of the official vcfpy docs, in docstrings, or even on the web in blog posts, articles, and such.

Snakemake uses `Sphinx <https://sphinx-doc.org>`_ for the user manual (that you are currently reading).
See :ref:`project_info-doc_guidelines` on how the documentation reStructuredText is used.


Submit Feedback
===============

The best way to send feedback is to file an issue at https://github.com/snakemake/snakemake/issues

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions are welcome :)

-----------------------
Pull Request Guidelines
-----------------------

To update the documentation, fix bugs or add new features you need to create a Pull Request
. A PR is a change you make to your local copy of the code for us to review and potentially integrate into the code base.

To create a Pull Request you need to do these steps:

1. Create a Github account.
2. Fork the repository.
3. Clone your fork locally.
4. Go to the created snakemake folder with :code:`cd snakemake`.
5. Create a new branch with :code:`git checkout -b <descriptive_branch_name>`.
6. Make your changes to the code or documentation.
7. Run :code:`git add .` to add all the changed files to the commit (to see what files will be added you can run :code:`git add . --dry-run`).
8. To commit the added files use :code:`git commit`. (This will open a command line editor to write a commit message. These should have a descriptive 80 line header, followed by an empty line, and then a description of what you did and why. To use your command line text editor of choice use (for example) :code:`export GIT_EDITOR=vim` before running :code:`git commit`).
9. Now you can push your changes to your Github copy of Snakemake by running :code:`git push origin <descriptive_branch_name>`.
10. If you now go to the webpage for your Github copy of Snakemake you should see a link in the sidebar called "Create Pull Request".
11. Now you need to choose your PR from the menu and click the "Create pull request" button. Be sure to change the pull request target branch to <descriptive_branch_name>!

If you want to create more pull requests, first run :code:`git checkout main` and then start at step 5. with a new branch name.

Feel free to ask questions about this if you want to contribute to Snakemake :)

------------------
Testing Guidelines
------------------

To ensure that you do not introduce bugs into Snakemake, you should test your code thoroughly.

Putting these tests into repeatable test cases ensures they can be checked on multiple platforms and Python versions.

Setup to run the test suite locally
===================================

Unit tests and regression tests are written to be run by `pytest <https://docs.pytest.org/en/stable/>`_.

Assuming you are on a Linux system, and have a working Conda installation, the standard way to run the tests is explained in the following.
With sufficient experience, it is however possible to run the tests in different setups, given the expected dependencies are installed.

1. Ensure your Conda version is at least 24.7.1 (Snakemake's minimum required Conda version).
   Check the output of ``conda --version`` and update conda if necessary.
2. Activate strict channel priorities (this is always a good idea when using Conda, see `here <https://conda-forge.org/docs/user/tipsandtricks/#using-multiple-channels>`__ for the rationale), by running ``conda config --set channel_priority strict``.
3. After checking out the branch you want to test, run these commands:

   .. code-block:: console

       $ conda env create -f test-environment.yml -n snakemake-testing
       $ conda activate snakemake-testing
       $ pip install -e .

   You may want to set a specific Python version by editing the constraint in ``test-environment.yml`` before doing this.
   Use of the ``-e``/``--editable`` option to ``pip`` will make your development version of Snakemake the one called when running Snakemake and all the unit tests. You only need to run the ``pip`` command once, not after each time you make code changes.

4. From the base Snakemake folder you may now run any specific test:

   .. code-block:: console

      $ pytest tests/tests.py::test_log_input

   You can also use the ``-k`` flag to select tests by substring match, rather than by the full name, and the ``--co`` option to preview which tests will be run. Try, for example:

   .. code-block:: console

      $ pytest --co tests/tests.py -k test_modules_all

Running the full test suite
===========================

If you simply run ``pytest`` in the top level directory it will scan for and attempt to run every test in the directory, but you will almost certainly get errors as not all tests are working and current.

The core test suite is run as a `GitHub action <https://docs.github.com/en/actions/about-github-actions/understanding-github-actions>`_ by the code under ``.github/workflows/main.yml``, so you should look in this file for the list of tests actually expected to pass in a regular test environment.

At the time of writing this text, the suite is:

.. code-block::

   tests/tests.py
   tests/tests_using_conda.py
   tests/test_expand.py
   tests/test_io.py
   tests/test_schema.py
   tests/test_linting.py
   tests/test_executor_test_suite.py
   tests/test_api.py
   tests/test_internals.py

Other tests in the directory may or may not work.

Warnings and oddities
=====================

You will likely see warnings related to deprecated functions in dependent libraries, especially botocore.

You may also get intermittent failures from tests that rely on external connectivity. The default test suite makes connections to multiple external services.

Tests that require singularity will be auto-skipped if no singularity or apptainer installation is available.
At the time of writing neither the ``singularity`` package on conda-forge nor the ``apptainer`` package are reliable, in that there are multiple failing tests on a standard Ubuntu system.
This is likely due to system security profiles that conda, being a non-root application, cannot change.
The Debian/Ubuntu ``singularity-container`` DEB package, which must be installed by the system administrator, does work.
The equivalent RPM package should also work on RedHat-type systems.

Tests in ``tests/test_api.py`` require a working ``git``.
This is not included in ``test-environment.yml`` as it's assumed you must have GIT installed to be working on the source code, but installing git into the conda environment should work if need be.

Depending on how the Snakemake code was downloaded and installed in the test environment, Snakemake may not be able to determine its own version and may think that it is version 0.
The existing unit tests should all cope with this, and in general you should avoid writing tests that rely on explicit version checks.

Continuous integration
======================

When creating a pull request for https://github.com/snakemake/snakemake, all tests will be automatically executed in a controlled environment by Github Actions.  


.. _project_info-doc_guidelines:

------------------------
Documentation Guidelines
------------------------

For the documentation, please adhere to the following guidelines:

- Put each sentence on its own line, this makes tracking changes through Git SCM easier.
- Provide hyperlink targets, at least for the first two section levels.
  For this, use the format ``<document_part>-<section_name>``, e.g., ``project_info-doc_guidelines``.
- Use the section structure from below.

::

    .. document_part-heading_1:

    =========
    Heading 1
    =========


    .. document_part-heading_2:

    ---------
    Heading 2
    ---------


    .. document_part-heading_3:

    Heading 3
    =========


    .. document_part-heading_4:

    Heading 4
    ---------


    .. document_part-heading_5:

    Heading 5
    ~~~~~~~~~


    .. document_part-heading_6:

    Heading 6
    :::::::::

.. _doc_setup:

-------------------
Documentation Setup
-------------------

For building the documentation, you have to install the Sphinx.
If you have already installed Conda, all you need to do is to create a
Snakemake development environment via

.. code-block:: console

    $ git clone git@github.com:snakemake/snakemake.git
    $ cd snakemake
    $ conda env create -f doc-environment.yml -n snakemake

You will also need to install your development version of Snakemake for the docs to be built correctly

.. code-block:: console

    $ pip install -e .

Then, the docs can be built with

.. code-block:: console

    $ conda activate snakemake
    $ cd docs
    $ make html
    $ make clean && make html  # force rebuild

Alternatively, you can use virtualenv.
The following assumes you have a working Python 3 setup.

.. code-block:: console

    $ git clone git@github.org:snakemake/snakemake.git
    $ cd snakemake/docs
    $ virtualenv -p python3 .venv
    $ source .venv/bin/activate
    $ pip install --upgrade -r requirements.txt

Afterwards, the docs can be built with

.. code-block:: console

    $ source .venv/bin/activate
    $ make html  # rebuild for changed files only
    $ make clean && make html  # force rebuild
