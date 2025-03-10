.. _snakefiles_configuration:

=============
Configuration
=============

Snakemake allows you to use configuration files for making your workflows more flexible and also for abstracting away direct dependencies to a fixed HPC cluster scheduler.


.. _snakefiles_standard_configuration:

----------------------
Standard Configuration
----------------------

Snakemake directly supports the configuration of your workflow.
A configuration is provided as a JSON or YAML file and can be loaded with:

.. code-block:: python

    configfile: "path/to/config.yaml"

The given path is interpreted relative to the working directory, not relative to the location of the snakefile that contains the statement.
The config file can be used to define a dictionary of configuration parameters and their values.
In case of YAML, the file can optionally be processed with `YTE <https://yte-template-engine.github.io>`_.
To activate this, you have to add the top-level key ``__use_yte__ = true`` to the YAML file.

In the workflow, the configuration is accessible via the global variable `config`, e.g.

.. code-block:: python

    rule all:
        input:
            expand("{sample}.{param}.output.pdf", sample=config["samples"], param=config["yourparam"])

If the `configfile` statement is not used, the config variable provides an empty dictionary.
In addition to the `configfile` statement, config values can be overwritten via the command line or the `snakemake.utils API <https://snakemake-api.readthedocs.io/en/stable/api_reference/snakemake_utils.html#snakemake.utils.update_config>`__, e.g.:

.. code-block:: console

    $ snakemake --config yourparam=1.5

Further, you can manually alter the config dictionary using any Python code **outside** of your rules. Changes made from within a rule won't be seen from other rules.
Finally, you can use the ``--configfile`` command line argument to overwrite values from the `configfile` statement.
Note that any values parsed into the ``config`` dictionary with any of above mechanisms are merged, i.e., all keys defined via a ``configfile``
statement, or the ``--configfile`` and ``--config`` command line arguments will end up in the final `config` dictionary, but if two methods define the same key, command line
overwrites the ``configfile`` statement.

For adding config placeholders into a shell command, Python string formatting syntax requires you to leave out the quotes around the key name, like so:

.. code-block:: python

    shell:
        "mycommand {config[foo]} ..."

.. _snakefiles_tabular_configuration:

---------------------
Tabular configuration
---------------------

It is usually advisable to complement YAML based configuration (see above) by a sheet based approach for meta-data that is of tabular form. For example, such
a sheet can contain per-sample information.
With the `Pandas library <https://pandas.pydata.org/>`_ such data can be read and used with minimal overhead, e.g.,

.. code-block:: python

    import pandas as pd

    samples = pd.read_table("samples.tsv").set_index("samples", drop=False)

reads in a table ``samples.tsv`` in TSV format and makes every record accessible by the sample name.
For details, see the `Pandas documentation <https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_table.html?highlight=read_table#pandas-read-table>`_.
A fully working real-world example containing both types of configuration can be found `here <https://github.com/snakemake-workflows/rna-seq-star-deseq2>`_.

---------------------
Environment variables
---------------------

Sometimes, it is not desirable to put configuration information into text files.
For example, this holds for secrets like access tokens or passwords.
Here, `environment variables <https://en.wikipedia.org/wiki/Environment_variable>`_ are the method of choice.
Snakemake allows to assert the existence of environment variables by adding a statement like:

.. code-block:: python

    envvars:
        "SOME_VARIABLE",
        "SOME_OTHER_VARIABLE"

When executing, Snakemake will fail with a reasonable error message if the variables ``SOME_VARIABLE`` and ``SOME_OTHER_VARIABLE`` are undefined.
Otherwise, it will take care of passing them to cluster and cloud environments. However, note that this does **not** mean that Snakemake makes them available e.g. in the jobs shell command.
Instead, for data provenance and reproducibility reasons, you are required to pass them explicitly to your job via the params directive, e.g. like this:

.. code-block:: python

    envvars:
        "SOME_VARIABLE"

    rule do_something:
        output:
             "test.txt"
        params:
            x=os.environ["SOME_VARIABLE"]
        shell:
            "echo {params.x} > {output}"


.. _snakefiles_config_validation:

----------
Validation
----------

With Snakemake 5.1, it is possible to validate both types of configuration via `JSON schemas <https://json-schema.org>`_.
The function ``snakemake.utils.validate`` takes a loaded configuration (a config dictionary or a Pandas data frame) and validates it with a given JSON schema.
Thereby, the schema can be provided in JSON or YAML format. Also, by using the defaults property it is possible to populate entries with default values. See `jsonschema FAQ on setting default values <https://python-jsonschema.readthedocs.io/en/latest/faq/>`_ for details.
In case of the data frame, the schema should model the record that is expected in each row of the data frame.
In the following example,

.. code-block:: python

  import pandas as pd
  from snakemake.utils import validate

  configfile: "config.yaml"
  validate(config, "config.schema.yaml")

  samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
  validate(samples, "samples.schema.yaml")


  rule all:
      input:
          expand("test.{sample}.txt", sample=samples.index)


  rule a:
      output:
          "test.{sample}.txt"
      shell:
          "touch {output}"

the schema for validating the samples data frame looks like this:

.. code-block:: yaml

  $schema: "https://json-schema.org/draft-06/schema#"
  description: an entry in the sample sheet
  properties:
    sample:
      type: string
      description: sample name/identifier
    condition:
      type: string
      description: sample condition that will be compared during differential expression analysis (e.g. a treatment, a tissue time, a disease)
    case:
      type: boolean
      default: true
      description: boolean that indicates if sample is case or control

  required:
    - sample
    - condition

Here, in case the case column is missing, the validate function will
populate it with True for all entries.

.. _snakefiles-peps:

-------------------------------------------
Configuring scientific experiments via PEPs
-------------------------------------------

Often scientific experiments consist of a set of samples (with optional subsamples), for which raw data and metainformation is known.
Instead of writing custom sample sheets as shown above, Snakemake allows to use `portable encapsulated project (PEP) <http://pep.databio.org>`_ definitions to configure a workflow.
This is done via a special directive `pepfile`, that can optionally complemented by a schema for validation (which is recommended for production workflows):

.. code-block:: python

    pepfile: "pep/config.yaml"
    pepschema: "schemas/pep.yaml"

    rule all:
        input:
            expand("{sample}.txt", sample=pep.sample_table["sample_name"])

    rule a:
        output:
            "{sample}.txt"
        shell:
            "touch {output}"

Using the ``pepfile`` directive leads to parsing of the provided PEP with `peppy <http://peppy.databio.org>`_.
The resulting project object is made globally available under the name ``pep``.
Here, we use it to aggregate over the set of sample names that is defined in the corresponding PEP.

**Importantly**, note that PEPs are meant to contain sample metadata and any global information about a project or experiment. 
They should **not** be used to encode workflow specific configuration options.
For those, one should always complement the pepfile with an ordinary :ref:`config file <snakefiles_standard_configuration>`.
The rationale is that PEPs should be portable between different data analysis workflows (that could be applied to the same data) and even between workflow management systems.
In other words, a PEP should describe everything needed about the data, while a workflow and its configuration should describe everything needed about the analysis that is applied to it.

^^^^^^^^^^^^^^^
Validating PEPs
^^^^^^^^^^^^^^^

Using the ``pepschema`` directive leads to an automatic parsing of the provided schema *and* PEP validation with the PEP validation tool -- `eido <http://eido.databio.org>`_. Eido schemas extend `JSON Schema <https://json-schema.org>`_ vocabulary to accommodate the powerful PEP features. Follow the `How to write a PEP schema <http://eido.databio.org/en/latest/writing-a-schema>`_ guide to learn more.

---------------------------
Configure Working Directory
---------------------------

All paths in the snakefile are interpreted relative to the directory snakemake is executed in. This behaviour can be overridden by specifying a workdir in the snakefile:

.. code-block:: python

    workdir: "path/to/workdir"

Usually, it is preferred to only set the working directory via the command line, because above directive limits the portability of Snakemake workflows.


.. _snakefiles-cluster_configuration:

---------------------------------------------
Cluster Configuration (not supported anymore)
---------------------------------------------

The previously supported cluster configuration has been replaced by configuration profiles (see :ref:`profiles`).
