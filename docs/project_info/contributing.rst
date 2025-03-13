.. _project_info-contributing:

************
Contributing
************

Contributions are welcome, and they are greatly appreciated!
Every little bit helps, and credit will always be given.
A detailed description of the Snakemake codebase and architecture can be found :ref:`here <codebase_intro>`.

You can contribute in many ways:

Types of contributions
======================


Report bugs or suggest enhancements
-----------------------------------

Report bugs or suggest enhancements at https://github.com/snakemake/snakemake/issues

If you are reporting a bug, follow the template and fill out the requested information, including:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug and ideally a test case.

If you are proposing am enhancement or a new feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions are welcome :).


Fix Bugs
--------

Look through the Github issues for bugs.
If you want to start working on a bug then please write short message on the issue tracker to prevent duplicate work.


Implement Features
------------------

Look through the Github issues for feature/enhancement requests.
If you want to start working on an issue then please write short message on the issue tracker to prevent duplicate work.

Contributing a plugin
---------------------

Currently, Snakemake supports executor plugins, storage plugins, and report plugins.
The `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog>`_ shows which plugins are available and how to contribute new ones.

Write Documentation
-------------------

Snakemake could always use more documentation, whether as part of the official docs, in docstrings, or even on the web in blog posts, articles, and such.

Snakemake uses `Sphinx <https://sphinx-doc.org>`_ for the user manual (that you are currently reading).
See :ref:`project_info-doc_guidelines` on how the documentation reStructuredText is used.



Pull Request Guidelines
=======================

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

.. _pixi-getting_started:

Development with ``pixi``
=========================

`pixi <https://pixi.sh/>`_ is a tool that is designed to help you manage 
your development environment.It acts as a drop-in replacement for
`conda <https://docs.conda.io/en/latest/>`_, offering:

- **Easy installation & Updating**: `install pixi <https://pixi.sh/latest/#installation>`_ 
  through many methods and for different shells.
  Updating ``pixi`` is as simple as ``pixi self-update``

- **Ease of Use**: A streamlined CLI (similar to Yarn or Cargo) for quick
  environment creation and management. Try commands like ``pixi init``,
  ``pixi add <package>``, or ``pixi run`` to see how intuitive it is.

- **Multiple Environments**: Define and switch between multiple sets of
  dependencies under one project.
  Pixi uses a ``feature`` system to compose an ``environment``.
  Think of ``features`` as a way to group dependencies and settings together.
  and an environment is a collection of features.
  This allows easy management of different environments for multiple use.
  See the ``pyproject.toml`` file for an example of how the ``test`` feature
  is used to define the ``dev``, ``py311`` and ``py312`` environments.

- **Cross-Platform Solving**: Target Linux, macOS, and Windows from a single
  config. Pixi resolves the correct packages for each platform and captures
  them in a lockfile for reproducible setups—no Docker needed.

- **Speed & Conda Compatibility**: Written in Rust, Pixi downloads and solves
  packages in parallel for faster operations. It uses the Conda ecosystem
  and channels, so you get the same packages with improved performance. In
  many cases, Pixi can outperform both Conda and Mamba.

To learn more, visit the `Pixi docs <https://pixi.sh>`__ or check out helpful
guides on `prefix.dev <https://prefix.dev/>`__. 

Testing Guidelines
==================

To ensure that you do not introduce bugs into Snakemake, you should test your code thoroughly.

Putting these tests into repeatable test cases ensures they can be checked on multiple platforms and Python versions.

Continuous integration
----------------------
For any pull request, all tests are automatically executed within Github Actions, providing feedback to you and the official development team whether the proposed changes are working as expected and do not hamper other functionality Snakemake provides.
However, it is useful to be able to run the tests locally, thereby being able to quickly debug any occurring failures.

Setup to run the test suite locally
-----------------------------------

Unit tests and regression tests are written to be run by `pytest <https://docs.pytest.org/en/stable/>`_.


.. _pixi-test-guide:

Testing Guide using ``pixi``
=============================

**Prerequisites**: Make sure you have ``pixi`` installed: See :ref:`pixi-getting_started`.

**Activate your environment**:
--------------------------------

There are a few environments you can use to run the tests.

The ``dev`` environment is most useful for overall development.
This environment will also install the ``docs`` and ``style`` features
which will allow you to also build documentation and run ``black``.

.. code-block:: console

    $ pixi shell -e dev

The ``py311`` and ``py312`` environments are what are used in the 
CI tests which isolate the Python version and the `test` dependencies.
Use this if you want to test your code against the same environment
as any failing CI tests.

.. code-block:: console

    $ pixi shell -e py311

**Run a comprehensive, simple, or single test**:
The test suite defines two types of tests via ``pixi tasks`` that you can run:

**test-all**: This task runs the comprehensive test suite, which includes 
*most* of the tests in the ``tests/`` directory.

.. code-block:: console

    $ pixi run test-all

**test-simple**: This task runs the main tests located in ``tests/tests.py``.

.. code-block:: console

    $ pixi run test-simple

**Single test**: You can also run a single test by using ``pytest`` 
directly with the test file and the test name.

.. code-block:: console

    $ pixi run pytest tests/tests.py::test_log_input

.. tip::
    This test suite is quite long, and can be run in parts similar to the 
    CI/CD tests which run it in 1/10 parts.

    To do so, you can use the ``--splits`` and ``--group`` flags to run
    a subset of the tests. For example, to run the first group of tests
    in a 10 part split:

    .. code-block:: console

        $ pixi run test-simple \
            --splits 10 \
            --group 1 \
            --splitting-algorithm=least_duration


Warnings and oddities
---------------------

You will likely see warnings related to deprecated functions in dependent libraries, especially botocore.

You may also get intermittent failures from tests that rely on external connectivity. The default test suite makes connections to multiple external services.

Tests that require singularity will be auto-skipped if no singularity or apptainer installation is available.
At the time of writing neither the ``singularity`` package on conda-forge nor the ``apptainer`` package are reliable, in that there are multiple failing tests on a standard Ubuntu system.
This is likely due to system security profiles that conda, being a non-root application, cannot change.
The Debian/Ubuntu ``singularity-container`` DEB package, which must be installed by the system administrator, does work.
The equivalent RPM package should also work on RedHat-type systems.

Depending on how the Snakemake code was downloaded and installed in the test environment, Snakemake may not be able to determine its own version and may think that it is version 0.
The existing unit tests should all cope with this, and in general you should avoid writing tests that rely on explicit version checks.


.. _project_info-doc_guidelines:

Documentation Guidelines
========================

For the documentation, please adhere to the following guidelines:

- Put each sentence on its own line, this makes tracking changes through Git SCM easier.
- Provide hyperlink targets, at least for the first two section levels.
  For this, use the format ``<document_part>-<section_name>``, e.g., ``project_info-doc_guidelines``.
- Use the `section structure recommended by Sphinx <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#sections>`_, which references the `recommendations in the Python Developer's Guide <https://devguide.python.org/documentation/markup/#sections>`_.
  Namely, the levels are:

::

    .. document_part-section_heading:

    ===============
    Section heading
    ===============


    .. document_part-subsection_heading:

    Subsection heading
    ------------------

    .. document_part-subsubsection_heading:

    Subsubsection heading
    ^^^^^^^^^^^^^^^^^^^^^

    .. document_part-paragraph_heading:

    Paragraph heading
    """""""""""""""""

.. _doc_setup:

Documentation Setup
-------------------

To get started, make sure you have ``pixi`` installed: 
See :ref:`pixi-getting_started`.
We use ``pixi`` to manage the docs environment and tasks to streamline
the developer experience.

.. code-block:: console

    $ ➜ pixi task list --environment docs
    Tasks that can run on this machine:
    -----------------------------------
    build-apidocs, build-docs, docs

    - build-apidocs   Build the API documentation in the apidocs/ directory
    - build-docs      Build the documentation in the docs/ directory
    - docs            Serve the documentation on http://localhost:8000 with live reload

**Test if the docs build**:
To only build the documentation, you can use the ``build-docs`` task.

.. code-block:: console

    $ pixi run build-docs

**Live server with auto-reload**:
To serve the documentation on a local server with live reload, 
use the ``docs`` task.

.. code-block:: console

    $ pixi run docs
    [sphinx-autobuild] Starting initial build
    [sphinx-autobuild] > python -m sphinx build docs/ docs/_build/html
    ...
    The HTML pages are in docs/_build/html.
    [sphinx-autobuild] Serving on http://0.0.0.0:8000
    [sphinx-autobuild] Waiting to detect changes...