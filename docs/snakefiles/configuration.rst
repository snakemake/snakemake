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

    configfile: "path/to/config.json"

The config file can be used to define a dictionary of configuration parameters and their values.
In the workflow, the configuration is accessible via the global variable `config`, e.g.

.. code-block:: python

    rule all:
        input:
            expand("{sample}.{param}.output.pdf", sample=config["samples"], param=config["yourparam"])

If the `configfile` statement is not used, the config variable provides an empty array.
In addition to the `configfile` statement, config values can be overwritten via the command line or the :ref:`api_reference_snakemake`, e.g.:

.. code-block:: console

    $ snakemake --config yourparam=1.5

Further, you can manually alter the config dictionary using any Python code **outside** of your rules. Changes made from within a rule won't be seen from other rules.
Finally, you can use the `--configfile` command line argument to overwrite values from the `configfile` statement.
Note that any values parsed into the `config` dictionary with any of above mechanisms are merged, i.e., all keys defined via a `configfile`
statement, or the `--configfile` and `--config` command line arguments will end up in the final `config` dictionary, but if two methods define the same key, command line
overwrites the `configfile` statement.

For adding config placeholders into a shell command, Python string formatting syntax requires you to leave out the quotes around the key name, like so:

.. code-block:: python

    shell:
        "mycommand {config[foo]} ..."

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

.. _snakefiles-cluster_configuration:

----------------------------------
Cluster Configuration (deprecated)
----------------------------------

While still being possible, **cluster configuration has been deprecated** by the introduction of :ref:`profiles`.

Snakemake supports a separate configuration file for execution on a cluster.
A cluster config file allows you to specify cluster submission parameters outside the Snakefile.
The cluster config is a JSON- or YAML-formatted file that contains objects that match names of rules in the Snakefile.
The parameters in the cluster config are then accessed by the ``cluster.*`` wildcard when you are submitting jobs.
Note that a workflow shall never depend on a cluster configuration, because this would limit its portability.
Therefore, it is also not intended to access the cluster configuration from **within** the workflow.

For example, say that you have the following Snakefile:

.. code-block:: python

    rule all:
        input: "input1.txt", "input2.txt"

    rule compute1:
        output: "input1.txt"
        shell: "touch input1.txt"

    rule compute2:
        output: "input2.txt"
        shell: "touch input2.txt"

This Snakefile can then be configured by a corresponding cluster config, say "cluster.json":


.. code-block:: json

    {
        "__default__" :
        {
            "account" : "my account",
            "time" : "00:15:00",
            "n" : 1,
            "partition" : "core"
        },
        "compute1" :
        {
            "time" : "00:20:00"
        }
    }

Any string in the cluster configuration can be formatted in the same way as shell commands, e.g. ``{rule}.{wildcards.sample}`` is formatted to ``a.xy`` if the rulename is ``a`` and the wildcard value is ``xy``.
Here ``__default__`` is a special object that specifies default parameters, these will be inherited by the other configuration objects. The ``compute1`` object here changes the ``time`` parameter, but keeps the other parameters from ``__default__``. The rule ``compute2`` does not have any configuration, and will therefore use the default configuration. You can then run the Snakefile with the following command on a SLURM system.

.. code-block:: console

    $ snakemake -j 999 --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n}  -t {cluster.time}"


For cluster systems using LSF/BSUB, a cluster config may look like this:

.. code-block:: json

    {
        "__default__" :
        {
            "queue"     : "medium_priority",
            "nCPUs"     : "16",
            "memory"    : 20000,
            "resources" : "\"select[mem>20000] rusage[mem=20000] span[hosts=1]\"",
            "name"      : "JOBNAME.{rule}.{wildcards}",
            "output"    : "logs/cluster/{rule}.{wildcards}.out",
            "error"     : "logs/cluster/{rule}.{wildcards}.err"
        },


        "trimming_PE" :
        {
            "memory"    : 30000,
            "resources" : "\"select[mem>30000] rusage[mem=30000] span[hosts=1]\"",
        }
    }

The advantage of this setup is that it is already pretty general by exploiting the wildcard possibilities that Snakemake provides via ``{rule}`` and ``{wildcards}``. So job names, output and error files all have reasonable and trackable default names, only the directies (``logs/cluster``) and job names (``JOBNAME``) have to adjusted accordingly.
If a rule named ``bamCoverage`` is executed with the wildcard ``basename = sample1``, for example, the output and error files will be ``bamCoverage.basename=sample1.out`` and ``bamCoverage.basename=sample1.err``, respectively.


---------------------------
Configure Working Directory
---------------------------

All paths in the snakefile are interpreted relative to the directory snakemake is executed in. This behaviour can be overridden by specifying a workdir in the snakefile:

.. code-block:: python

    workdir: "path/to/workdir"

Usually, it is preferred to only set the working directory via the command line, because above directive limits the portability of Snakemake workflows.
