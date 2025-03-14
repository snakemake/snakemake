
.. _codebase_intro:

Codebase and architecture
=========================

Snakemake is organized as a combination of a main package, a set of plugin interface packages for various functionalities (e.g. execution backends, storage, reporting), and plugin packages implementing those interfaces.
The development of the Snakemake main package as well as the plugin interfaces is hosted by the Snakemake GitHub organization (https://github.com/snakemake) and maintained by the :ref:`Snakemake core team <maintainers>`.
Many plugins are hosted by and developed within the Snakemake GitHub organization as well.
However, both Snakemake's plugin detection and the `plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog/>`__ are designed to work seamlessly with plugins developed outside of the Snakemake GitHub organization and such contributions are highly encouraged.
Visually, Snakemake's architecture can be summarized as follows (note that some plugin interfaces are still in development and not yet available to the public):

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
The ``ioutils`` module (`snakemake/ioutils <https://github.com/snakemake/snakemake/blob/main/snakemake/ioutils>`__) implements semantic helper functions for handling input and output files as well as non-file parameters in the workflow.

linting
"""""""
The ``linting`` module (`snakemake/linting <https://github.com/snakemake/snakemake/blob/main/snakemake/linting>`__) implements static code analysis functionality for giving hints on discovered anti-patterns in a given workflow definition.

script
""""""
The ``script`` module (`snakemake/script <https://github.com/snakemake/snakemake/blob/main/snakemake/script>`__) implements Snakemake's integration with scripting languages like R or Julia.


Plugins
^^^^^^^

Various functionalities of Snakemake are organized via plugins (e.g. supporting different execution backends or remote storage).
For each type of plugin, Snakemake offers an interface package.
Interface packages provide a stable API for the communication between the main Snakemake and the plugin.
They are independently (and semantically) versioned and aim to avoid breaking changes, so that plugins will rarely have to be modified and should remain compatible with a wide range of Snakemake versions for a long time.

Interfaces
""""""""""

Interface packages offer a set of strongly typed abstract base classes that have to be implemented by the respective plugins, thereby rigorously defining what functionality a plugin has to offer and to which structures of the main Snakemake it has access.
Each interface package comes with a test suite, that tests any implementing plugin in terms of (ideally) all kinds of plugin usages.

Naming
""""""

Plugins have to follow a naming convention, controlling how they are findable via `pypi <https://pypi.io>`__.
The convention state that the plugin has to be named as ``snakemake-<type>-plugin-<name>`` with ``<type>`` being the type of plugin (e.g. ``storage``) and ``<name>`` being the name of the plugin (e.g. ``s3``).

Plugin catalog
""""""""""""""

Leveraging this standardized naming, Snakemake automatically collects all plugins from pypi and presents them in the `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog>`__.
The catalog automatically generates usage documentation for each plugin, gives credit to authors, and provides links to the respective github repositories.
Plugin authors can extend the catalog documentation for their plugin by providing markdown files ``docs/intro.md`` and ``docs/further.md`` in the plugin repository in order to show an introductory paragraph as well as extensive non-standard documentation (like usage examples or other plugin specific information) in the catalog.

Scaffolding
"""""""""""

Via `poetry <https://github.com/snakemake/poetry-snakemake-plugin>`__, plugins can be automatically scaffolded, leading to all files for plugin implementation, testing and package building being generated as skeletons.
The developer then only needs to implement extensively annotated abstract methods of base classes provided by the respective interface package.
After pushing the plugin code into a new GitHub repository, testing, release-automation, and pypi upload then work out of the box.

Assuming that the plugin type to create is given as bash variable ``$type`` below, given that poetry is available as a command, the following procedure should be followed for scaffolding a new plugin:

.. code-block:: bash

    # Install latest version of the snakemake poetry plugin
    poetry self add poetry-snakemake-plugin@latest

    # Create a new poetry project via
    poetry new snakemake-$type-plugin-myfancyplugin

    cd snakemake-$type-plugin-myfancyplugin

    poetry scaffold-snakemake-$type-plugin

    # Next, edit the scaffolded code according to your needs, and publish
    # the resulting plugin into a github repository. The scaffold command also 
    # creates github actions workflows that will immediately start to check and test
    # the plugin.

An example class (in this case for storage plugins) created by the scaffold command would be the following:

.. code-block:: python

    # Required:
    # Implementation of your storage provider
    # This class can be empty as the one below.
    # You can however use it to store global information or maintain e.g. a connection
    # pool.
    class StorageProvider(StorageProviderBase):
        # For compatibility with future changes, you should not overwrite the __init__
        # method. Instead, use __post_init__ to set additional attributes and initialize
        # further stuff.

        def __post_init__(self):
            # This is optional and can be removed if not needed.
            # Alternatively, you can e.g. prepare a connection to your storage backend here.
            # and set additional attributes.
            pass

        @classmethod
        def example_queries(cls) -> List[ExampleQuery]:
            # Return example queries with description for this storage provider (at
            # least one).
            ...

        def rate_limiter_key(self, query: str, operation: Operation) -> Any:
            # Return a key for identifying a rate limiter given a query and an operation.

            # This is used to identify a rate limiter for the query.
            # E.g. for a storage provider like http that would be the host name.
            # For s3 it might be just the endpoint URL.
            ...

        def default_max_requests_per_second(self) -> float:
            # Return the default maximum number of requests per second for this storage
            # provider.
            ...

        def use_rate_limiter(self) -> bool:
            # Return False if no rate limiting is needed for this provider.
            ...

        @classmethod
        def is_valid_query(cls, query: str) -> StorageQueryValidationResult:
            # Return whether the given query is valid for this storage provider.
            # Ensure that also queries containing wildcards (e.g. {sample}) are accepted
            # and considered valid. The wildcards will be resolved before the storage
            # object is actually used.
            ...

Once all methods of all scaffolded classes are implemented, the plugin is ready to be tested.

Continuous testing is conducted via Github Actions, defined in the file ``.github/workflows/ci.yml``.
In case the testing needs additional software or services to be deployed for the plugin to be tested, this can happen inside that file, prior to the step that invokes pytest.
