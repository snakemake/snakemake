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

To have integration tests run automatically when committing code changes to Github, you need to sign up on wercker.com and register a user.

The easiest way to run your development version of Snakemake is perhaps to go to the folder containing your local copy of Snakemake and call:

.. code-block:: console

    $ conda env create -f test-environment.yml -n snakemake-testing
    $ conda activate snakemake-testing
    $ pip install -e .

This will make your development version of Snakemake the one called when running snakemake. You do not need to run this command after each time you make code changes.

From the base snakemake folder you call :code:`pytest` to run all the tests, or choose one specific test:

.. code-block:: console

   $ pytest
   $ pytest tests/tests.py::test_log_input

If you introduce a new feature you should add a new test to the tests directory. See the folder for examples.

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
