.. _project_info-contributing:

************
Contributing
************

Contributions are welcome, and they are greatly appreciated!
Every little bit helps, and credit will always be given.

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
---------------------------

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
---------------------

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


.. _project_info-doc_guidelines:

Documentation Guidelines
========================

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

.. _codebase_intro:

The Snakemake codebase
======================

Snakemake is organized as a combination of a main package, a set of plugin interface packages for various functionalities (e.g. execution backends, storage, reporting), and plugin packages implementing those interfaces.
The development of the Snakemake main package as well as the plugin interfaces is hosted by the Snakemake GitHub organization (https://github.com/snakemake) and maintained by the :ref:`Snakemake core team <maintainers>`.
Many plugins are hosted by and developed within the Snakemake GitHub organization as well.
However, both Snakemake's plugin detection and the `plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog/>`__ are designed to work seamlessly with plugins developed outside of the Snakemake GitHub organization and such contributions are highly encouraged.
Visually, Snakemake's architecture can be summarized as follows:

.. image:: img/architecture.svg
    :alt: Snakemake architecture

The main package
----------------

Snakemake's main package can be partitioned into three levels.

#. The user-facing level, consisting of an API and a command line interface that uses the API under the hood.
#. The language level, implementing the language parser and thereby defining the workflow definition syntax.
#. The core level, implementing Snakemake's interpretation of rules defined via the language, the inference of jobs and their dependencies, as well as the actual execution of workflows as a combined usage of the various plugin types.

The user-facing level
^^^^^^^^^^^^^^^^^^^^^

Most users will interact with the command line interface of Snakemake.
The command line interface of Snakemake is defined in `snakemake/cli.py <https://github.com/snakemake/snakemake/blob/main/snakemake/cli.py>`__ using Python's argparse module.
The module is responsible for parsing command-line arguments, setting up logging, and invoking the appropriate Snakemake API calls based on the provided arguments. It defines various argument groups for execution, grouping, reports, notebooks, utilities, output, behavior, remote execution, software deployment, conda, apptainer/singularity, environment modules, and internal use.

The main steps in `cli.py` include:
#. **Argument Parsing**: Using argparse to define and parse command-line arguments.
#. **Profile Handling**: Managing profiles for setting default command-line arguments.
#. **Logging Setup**: Configuring log handlers based on the provided arguments.
#. **API Invocation**: Converting parsed arguments into API calls to execute the workflow, generate reports, or perform other tasks.

The API (application programming interface) of Snakemake is defined in `snakemake/api.py <https://github.com/snakemake/snakemake/blob/main/snakemake/api.py>`__.
It provides a programmatic interface to Snakemake, allowing users to interact with Snakemake workflows through Python code.
It defines classes and methods for setting up and executing workflows, managing resources, and handling various settings.

The main components include:
#. ``SnakemakeApi``: The main entry point for interacting with Snakemake workflows programmatically. It provides methods for setting up workflows, deploying sources, and handling exceptions.
   Via the method ``SnakemakeApi.workflow()`` a workflow can be instantiated from a given snakefile, thereby returning a ``WorkflowApi`` object).
#. ``WorkflowApi``: A class for managing workflow-specific settings and operations, such as creating DAGs, linting workflows, and executing workflows.
   Via the method ``WorkflowApi.dag()`` a Directed Acyclic Graph (DAG) can be instantiated from the workflow, thereby returning a ``DAGApi`` object).
#. ``DAGApi``: A class for managing Directed Acyclic Graph (DAG) settings and executing workflows with specific executors and settings.

All parts of the API are designed to be used in combination with a set of `data classes <https://docs.python.org/3/library/dataclasses.html>`__ and `enumerations <https://docs.python.org/3/library/enum.html>`__ that define various settings in a type safe and readable way.
These classes are in part found under `snakemake/settings <https://github.com/snakemake/snakemake/blob/main/snakemake/settings>`__, and in part (where necessary) defined in the plugin interface packages.

The language level
^^^^^^^^^^^^^^^^^^

Snakemake offers a domain specific language (DSL) for defining workflows.
The syntax of the DSL is primarily defined via the implementation of the parser in `snakemake/parser.py <https://github.com/snakemake/snakemake/blob/main/snakemake/parser.py>`__.
A key feature of Snakemake's power is the fact that the language extends the syntax of Python with *directives* to define rules and other workflow specific controls.
Technically, this is implemented as a hierarchy of `Mealy machines <https://en.wikipedia.org/wiki/Mealy_machine>`__, each of which is responsible for one of the directives offered by the Snakemake DSL.
The parser first tokenizes a given snakefile using Python's builtin tokenizer, and then translates the tokens (that contain DSL tokens) into plain python tokens which translate the DSL specification into a plain Python specification of the workflow and all its items (e.g. the definition of rules).
The latter are then fed into the Python interpreter, thereby building up the workflow with all the defined rules.
During the token translation process, the Mealy machine hierarchy starts in a state accepting (and outputting) any Python token, and switches to appropriate sub-machines whenever a DSL-specifiy keyword is encountered.
The Mealy machines are implemented as Python classes, using abstract base classes to define common functionality.
The generated Python code invokes methods of the ``snakemake.workflow.Workflow`` (see :ref:`codebase_core`) class to build up the workflow.

.. _codebase_core:

The core level
^^^^^^^^^^^^^^

The core level of Snakemake is responsible for interpreting the rules defined via the language, inferring jobs and their dependencies, and executing workflows. The most important components of the core level include:

Workflow
""""""""
The ``Workflow`` class (`snakemake/workflow.py <https://github.com/snakemake/snakemake/blob/main/snakemake/workflow.py>`__) is the central class representing a Snakemake workflow. It manages the rules, config, resources, and execution settings. The `Workflow` class is responsible for parsing the workflow definition, creating the Directed Acyclic Graph (DAG) of jobs, and orchestrating the execution of the workflow. It interacts with other core components like `Rule`, `DAG`, and `Scheduler` to manage the workflow execution.

Rule
""""
The `Rule` class (`snakemake/rules.py <https://github.com/snakemake/snakemake/blob/main/snakemake/rules.py>`__) represents a single rule in the workflow. A `Rule` defines the input, output, and parameters for a specific step in the workflow. It also includes directives for resources, conda environments, and containerization. The `Rule` class interacts with the `Workflow` class to define the workflow structure and with the `DAG` class to manage job dependencies.

DAG
"""
The `DAG` class (`snakemake/dag.py <https://github.com/snakemake/snakemake/blob/main/snakemake/dag.py>`__) represents the Directed Acyclic Graph (DAG) of jobs in the workflow. It is responsible for inferring the order of job execution, detecting cycles, and managing job dependencies. The `DAG` class interacts with the `Workflow`, `Rule`, and `Scheduler` classes to manage the workflow execution.

Persistence
"""""""""""
The `Persistence` class (`snakemake/persistence.py <https://github.com/snakemake/snakemake/blob/main/snakemake/persistence.py>`__) manages the persistent storage of metadata and provenance information. The `Persistence` class interacts with the `Workflow` and `DAG` classes to manage workflow state.

Scheduler
"""""""""
The `Scheduler` class (`snakemake/scheduler.py <https://github.com/snakemake/snakemake/blob/main/snakemake/scheduler.py>`__) is responsible for scheduling jobs for execution. It uses the DAG to determine the order of job execution, taking into account resource constraints and job priorities. The `Scheduler` class interacts with the `Workflow`, `DAG`, and `Rule` classes to manage job scheduling.

PathModifier
""""""""""""
The `PathModifier` class (`snakemake/path_modifier.py <https://github.com/snakemake/snakemake/blob/main/snakemake/path_modifier.py>`__) is a utility class for handling path modifications such as the handling of remote storage or module imports. It ensures that paths are correctly managed and modified according to the workflow's requirements. The `PathModifier` class interacts with the `Workflow` and `_IOFile` classes to manage file paths.

Sourcecache
"""""""""""
The `Sourcecache` class (`snakemake/sourcecache.py <https://github.com/snakemake/snakemake/blob/main/snakemake/sourcecache.py>`__) handles the caching of source files. It ensures that remote source files are efficiently managed and reused across workflow executions. The `Sourcecache` class interacts with the `Workflow` and `_IOFile` classes to manage source files.

Besides these central classes, the following modules add additional functionality:

ioutils
"""""""
The ``ioutils`` module (`snakemake/ioutils <https://github.com/snakemake/snakemake/blob/main/snakemake/ioutils>`__) implements semantic helper functions functions for handling input and output files as well as non-file parameters in the workflow.
