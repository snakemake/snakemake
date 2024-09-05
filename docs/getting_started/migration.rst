.. _migration:

====================================
Migration between Snakemake versions
====================================

Snakemake is meant to remain backwards compatible as much as possible.
However, sometimes, very rarely, we remove old almost unused features that have since then
been replaced by new ones (so far, this happened only once, for Snakemake 8).
Sometimes, new features are added that do not require, but make it strongly advisable to adapt workflows (e.g. because the new features provide a better user or recipient experience).

Below are migration hints for particular Snakemake versions.

Migrating to Snakemake 8
------------------------

Workflow definitions
^^^^^^^^^^^^^^^^^^^^

Snakemake 8 removes the support for three syntactical elements, which are all officially deprecated since multiple major releases:

* Support for marking output files as ``dynamic`` has been removed. You should instead use :ref:`checkpoints <snakefiles-checkpoints>`.
* Support for the ``version`` directive has been removed. You should use the :ref:`conda <integrated_package_management>` or :ref:`container <apptainer>` integration instead.
* Support for the ``subworkflow`` directive has been removed. You should use the :ref:`module directive <snakefiles-modules>` instead, which provides the same functionality in a more general way.

In addition, we have moved the former remote provider functionality into so called :ref:`storage plugins <storage-support>`.
Most of the old remote providers have been migrated into the new storage plugins
(see the `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog>`__.).
Two former remote providers have been migrated into Snakemake wrappers instead, namely
the NCBI and EGA remote providers, which are now replaced by the
`entrez/efetch <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/entrez/efetch.html>`_ and
the `ega <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/ega/fetch.html>`_ wrappers.
As of writing, the Snakemake storage plugin for xrootd (see `here <https://github.com/snakemake/snakemake-storage-plugin-xrootd>`__) does not yet pass the CI tests. Any help would be greatly appreciated.


Command line interface
^^^^^^^^^^^^^^^^^^^^^^

The command line interface of Snakemake 8 has a lot of new options which are best explored using::

    snakemake --help

Moreover, some options have been renamed:

* All the execution backends have been moved into plugins. When you used e.g. ``--kubernetes`` and corresponding options before, you should now use ``--executor kubernetes`` and check the `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/kubernetes.html>`_ for the new options. The same holds for all other execution backends, see `here <https://snakemake.github.io/snakemake-plugin-catalog/index.html>`__.
* The ``--use-conda`` and ``--use-singularity`` options are deprecated. Instead you should now use ``--software-deployment-method conda`` or ``--software-deployment-method apptainer`` or ``--software-deployment-method conda apptainer`` if you need both.
* There is a new executor plugin for `Google Cloud Batch <https://cloud.google.com/batch/docs/get-started>`_.
  This is meant as a replacement for the old Google Life Sciences executor. 
  The new executor is called ``googlebatch`` and can be used with ``--executor googlebatch``. 
  Please check out the documentation of the plugin in the `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/googlebatch.html>`__. 
  Note that in principle it is fine to re-add google-lifesciences support as a plugin as well. 
  We even have skeleton code for this `here <https://github.com/snakemake/snakemake-executor-plugin-google-lifesciences>`__. 
  Any help with getting this tested and released despite the fact that google lifesciences will be shut down this year would still be valued.

.. list-table:: Interface comparison
   :widths: 15 30 15 30 10
   :header-rows: 1
   :align: center

   * - Interface in 7.32
     - Interface description in 7.32
     - Interface in 8.0.1
     - Interface description in 8.0.1
     - Change introduction
   * - preemptible
     -
     -
     -
     -
   * - --preemption-default PREEMPTION_DEFAULT
     -
                        A preemptible instance can be requested when using the
                        Google Life Sciences API. If you set a --preemption-
                        default, all rules will be subject to the default.
                        Specifically, **this integer is the number of restart
                        attempts** that will be made given that the instance is
                        killed unexpectedly. Note that preemptible instances
                        have a maximum running time of 24 hours. If you want
                        to set preemptible instances for only a subset of
                        rules, use --preemptible-rules instead. (default:
                        None)
     - --preemptible-retries PREEMPTIBLE_RETRIES
     -
                        **Number of retries** that shall be made in order to
                        finish a job from of rule that has been marked as
                        preemptible via the --preemptible-rules setting.
                        (default: None)
     - Renamed
   * - --preemptible-rules PREEMPTIBLE_RULES [PREEMPTIBLE_RULES ...]
     -
                        A preemptible instance can be requested when using the
                        Google Life Sciences API. If you want to use these
                        instances for a subset of your rules, you can use
                        --preemptible-rules and then specify a list of rule
                        and integer pairs, where each integer indicates the
                        number of restarts to use for the rule's instance in
                        the case that the instance is terminated unexpectedly.
                        --preemptible-rules can be used in combination with
                        --preemption-default, and will take priority. Note
                        that preemptible instances have a maximum running time
                        of 24. If you want to apply a consistent number of
                        retries across all your rules, use --preemption-
                        default instead. Example: snakemake --preemption-
                        default 10 --preemptible-rules map_reads=3
                        call_variants=0 (default: None)
     - --preemptible-rules [PREEMPTIBLE_RULES ...]
     -
                        Define which rules shall use a preemptible machine
                        which can be prematurely killed by e.g. a cloud
                        provider (also called spot instances). This is
                        currently only supported by the Google Life Sciences
                        executor and ignored by all other executors. If no
                        rule names are provided, all rules are considered to
                        be preemptible. The (default: None)
     - Renamed
   * - list-rules
     -
     -
     -
     -
   * - --list, -l
     -
                        Show available rules in given Snakefile. (default:
                        False)
     - **--list-rules**, --list, -l
     -
                        Show available rules in given Snakefile. (default:
                        False)
     - New alias: --list-rules
   * - list-changes
     -
     -
     -
     -
   * - --list-version-changes, --lv
     -
                        List all output files that have been created with a
                        different version (as determined by the version
                        keyword). (default: False)
     -
     -
     - Deprecated: It seems due to the deprecation of ``version`` directive
   * - --list-code-changes, --lc
     -
                        List all output files for which the rule body (run or
                        shell) have changed in the Snakefile. (default: False)
     - --list-changes {params,input,code}, --lc {params,input,code}
     -
                        List all output files for which the rule body (run or
                        shell) have changed in the Snakefile. (default: None)
     - Redesigned: Please use params such as ``--list-changes params,input,code`` instead of ``--list-code-changes``, ``--list-input-changes``, or ``--list-params-changes``
   * - bash-completion
     -
     -
     -
     -
   * - --bash-completion
     -
                        Output code to register bash completion for snakemake.
                        Put the following in your .bashrc (including the
                        accents): `snakemake --bash-completion` or issue it in
                        an open terminal session. (default: False)
     -
     -
     - Unsupported?
   * - deploy-sources
     -
     -
     -
     -
   * -
     -
     - --deploy-sources QUERY CHECKSUM
     -
                        Deploy sources archive from given storage provider
                        query to the current working sdirectory and control
                        for archive checksum to proceed. Meant for internal
                        use only. (default: None)
     -
   * - reason
     -
     -
     -
     -
   * - --reason, -r
     -
                        Print the reason for each executed rule (deprecated,
                        always true now). (default: False)
     -
     -
     - Deprecated: Drop it and don't worry about anything
   * - gui
     -
     -
     -
     -
   * - --gui [PORT]
     -
                        Serve an HTML based user interface to the given
                        network and port e.g. 168.129.10.15:8000. By default
                        Snakemake is only available in the local network
                        (default port: 8000). To make Snakemake listen to all
                        ip addresses add the special host address 0.0.0.0 to
                        the url (0.0.0.0:8000). This is important if Snakemake
                        is used in a virtualised environment like Docker. If
                        possible, a browser window is opened. (default: None)

     -
     -
     - Unsupported?
   * - stats
     -
     -
     -
     -
   * - --stats FILE
     -
                        Write stats about Snakefile execution in JSON format
                        to the given file. (default: None)
     -
     -
     - Unsupported?
   * - file storage
     -
     -
     -
     -
   * -
     -
     - --unneeded-temp-files FILE [FILE ...]
     -
                        Given files will not be uploaded to storage and
                        immediately deleted after job or group job completion.
                        (default: frozenset())
     -
   * - --keep-remote
     -
                        Keep local copies of remote input files. (default:
                        False)
     - --keep-storage-local-copies
     -
                        Keep local copies of remote input files. (default:
                        False)
     - Renamed
   * - --keep-target-files
     -
                        Do not adjust the paths of given target files relative
                        to the working directory. (default: False)
     -  --target-files-omit-workdir-adjustment
     -
                        Do not adjust the paths of given target files relative
                        to the working directory. (default: False)
     - Renamed
   * - seconds-between-status-checks
     -
     -
     -
     -
   * -
     -
     - --seconds-between-status-checks SECONDS_BETWEEN_STATUS_CHECKS
     -
                        Number of seconds to wait between two rounds of status
                        checks. (default: 10)
     -
   * - remote storage
     -
     -
     -
     -
   * - --default-remote-provider {S3,GS,FTP,SFTP,S3Mocked,gfal,gridftp,iRODS,AzBlob,XRootD}
     -
                        Specify default remote provider to be used for all
                        input and output files that don't yet specify one.
                        (default: None)
     - --default-storage-provider DEFAULT_STORAGE_PROVIDER
     -
                        Specify default storage provider to be used for all
                        input and output files that don't yet specify one
                        (e.g. 's3'). See
                        https://snakemake.github.io/snakemake-plugin-catalog
                        for available storage provider plugins. (default:
                        None)
     - Renamed:
                        See
                        https://snakemake.github.io/snakemake-plugin-catalog
                        for available storage provider plugins.
   * - --default-remote-prefix DEFAULT_REMOTE_PREFIX
     -
                        Specify prefix for default remote provider. E.g. a
                        bucket name. (default: )
     - --default-storage-prefix DEFAULT_STORAGE_PREFIX
     -
                        Specify prefix for default storage provider. E.g. a
                        bucket name. (default: )
     - Renamed
   * -
     -
     - --local-storage-prefix LOCAL_STORAGE_PREFIX
     -
                        Specify prefix for storing local copies of storage
                        files and folders. By default, this is a hidden
                        subfolder in the workdir. It can however be freely
                        chosen, e.g. in order to store those files on a local
                        scratch disk. (default: .snakemake/storage)
     -
   * - shared-fs
     -
     -
     -
     -
   * - --no-shared-fs
     -
                        Do not assume that jobs share a common file system.
                        When this flag is activated, Snakemake will assume
                        that the filesystem on a cluster node is not shared
                        with other nodes. For example, this will lead to
                        downloading remote files on each cluster node
                        separately. Further, it won't take special measures to
                        deal with filesystem latency issues. This option will
                        in most cases only make sense in combination with
                        --default-remote-provider. Further, when using
                        --cluster you will have to also provide --cluster-
                        status. Only activate this if you know what you are
                        doing. (default: False)
     - --shared-fs-usage {input-output,persistence,software-deployment,source-cache,sources,storage-local-copies,none} [{input-output,persistence,software-deployment,source-cache,sources,storage-local-copies,none} ...]
     -
                        Set assumptions on shared filesystem for non-local
                        workflow execution. To disable any sharing via the
                        filesystem, specify 'none'. Usually, the executor
                        plugin sets this to a correct default. However,
                        sometimes it is worth tuning this setting, e.g. for
                        optimizing cluster performance. For example, when
                        using '--default-storage-provider fs' together with a
                        cluster executor like slurm, you might want to set '--
                        shared-fs-usage persistence software-deployment
                        sources source-cache', such that software deployment
                        and data provenance will be handled by NFS but input
                        and output files will be handled exclusively by the
                        storage provider. (default:
                        frozenset({<SharedFSUsage.SOFTWARE_DEPLOYMENT: 2>,
                        <SharedFSUsage.INPUT_OUTPUT: 1>,
                        <SharedFSUsage.PERSISTENCE: 0>,
                        <SharedFSUsage.SOURCES: 3>,
                        <SharedFSUsage.SOURCE_CACHE: 5>,
                        <SharedFSUsage.STORAGE_LOCAL_COPIES: 4>}))
     - Redesigned: Please change ``--no-shared-fs`` to ``--shared-fs-usage none``
   * -
     -
     - --job-deploy-sources
     -
                        Whether the workflow sources shall be deployed before
                        a remote job is started. Only applies if --no-shared-
                        fs is set or executors are used that imply no shared
                        FS (e.g. the kubernetes executor). (default: False)
     - (Clearer description needed)
   * - greediness
     -
     -
     -
     -
   * - --greediness GREEDINESS
     -
                        Set the greediness of scheduling. This value between 0
                        and 1 determines how careful jobs are selected for
                        execution. The default value (1.0) provides the best
                        speed and still acceptable scheduling quality.
                        (default: None)
     - --scheduler-greediness SCHEDULER_GREEDINESS, --greediness SCHEDULER_GREEDINESS
     -
                        Set the greediness of scheduling. This value between 0
                        and 1 determines how careful jobs are selected for
                        execution. The default value (1.0) provides the best
                        speed and still acceptable scheduling quality.
                        (default: None)
     - Renamed
   * - debug
     -
     -
     -
     -
   * - --overwrite-shellcmd OVERWRITE_SHELLCMD
     -
                        Provide a shell command that shall be executed instead
                        of those given in the workflow. This is for debugging
                        purposes only. (default: None)
     -
     -
     - Deprecated
   * -  --mode {0,1,2}
     -
                        Set execution mode of Snakemake (internal use only).
                        (default: 0)
     - --mode {default,remote,subprocess}
     -
                        Set execution mode of Snakemake (internal use only).
                        (default: default)
     - Redesigned: use string instead of integer
   * - APPTAINER/SINGULARITY
     -
     -
     -
     -
   * - --use-singularity
     -
                        If defined in the rule, run job within a singularity
                        container. If this flag is not set, the singularity
                        directive is ignored. (default: False)
     - --use-apptainer, --use-singularity
     -
                        If defined in the rule, run job within a
                        apptainer/singularity container. If this flag is not
                        set, the singularity directive is ignored. (default:
                        False)
     - New alias (more general usage)
   * - --singularity-prefix DIR
     -
                        Specify a directory in which singularity images will
                        be stored. If not supplied, the value is set to the
                        '.snakemake' directory relative to the invocation
                        directory. If supplied, the ``--use-singularity`` flag
                        must also be set. The value may be given as a relative
                        path, which will be extrapolated to the invocation
                        directory, or as an absolute path. (default: None)
     - --apptainer-prefix DIR, --singularity-prefix DIR
     -
                        Specify a directory in which apptainer/singularity
                        images will be stored.If not supplied, the value is
                        set to the '.snakemake' directory relative to the
                        invocation directory. If supplied, the ``--use-
                        apptainer`` flag must also be set. The value may be
                        given as a relative path, which will be extrapolated
                        to the invocation directory, or as an absolute path.
                        (default: None)
     - New alias (more general usage)
   * - --singularity-args ARGS
     -
                        Pass additional args to singularity. (default: )
     - --apptainer-args ARGS, --singularity-args ARGS
     -
                        Pass additional args to apptainer/singularity.
                        (default: )
     - New alias (more general usage)
   * - --cleanup-containers
     -
                        Remove unused (singularity) containers (default:
                        False)
     - --container-cleanup-images
     -
                        Remove unused containers (default: False)
     - New alias (more general usage)
   * - precommand
     -
     -
     -
     -
   * - --precommand PRECOMMAND
     -
                        Any command to execute before snakemake command **on AWS
                        cloud** such as wget, git clone, unzip, etc. This is
                        used with --tibanna.Do not include input/output
                        download/upload commands - file transfer between S3
                        bucket and the run environment (container) is
                        automatically handled by Tibanna. (default: None)
     - --precommand PRECOMMAND
     -
                        Only used in case of remote execution. Command to be
                        executed before Snakemake executes each job on the
                        remote compute node. (default: None)
     - Redesigned: more general usage
   * - software-deployment-method
     -
     -
     -
     -
   * -
     -
     - --software-deployment-method {apptainer,conda,env-modules} [{apptainer,conda,env-modules} ...], --deployment-method {apptainer,conda,env-modules} [{apptainer,conda,env-modules} ...], --deployment {apptainer,conda,env-modules} [{apptainer,conda,env-modules} ...]
     -
                        Specify software environment deployment method.
                        (default: set())
     - New designed
   * - executor
     -
     -
     -
     -
   * - --cluster CMD, (may be --touch, --dryrun, ..., ?)
     -
     - --executor {cluster-generic,local,dryrun,touch}, -e {cluster-generic,local,dryrun,touch}
     -
                        Specify a custom executor, available via an executor
                        plugin: snakemake_executor_<name> (default: None)
     - New designed: Now if you want to use ``--cluster CMD``, please use ``--executor cluster-generic --cluster-generic-submit-cmd CMD`` instead.
        Note you should install ``cluster-generic`` using command ``pip install snakemake-executor-cluster-generic``
   * - --cluster CMD
     -
                        Execute snakemake rules with the given submit command,
                        e.g. qsub. Snakemake compiles jobs into scripts that
                        are submitted to the cluster with the given command,
                        once all input files for a particular job are present.
                        The submit command can be decorated to make it aware
                        of certain job properties (name, rulename, input,
                        output, params, wildcards, log, threads and
                        dependencies (see the argument below)), e.g.: $
                        snakemake --cluster 'qsub -pe threaded {threads}'.
                        (default: None)
     -  --cluster-generic-submit-cmd VALUE
        (Requires the cluster-generic plugin)
     -
                        Command for submitting jobs (default:
                        <dataclasses._MISSING_TYPE object at 0x7fc423088680>)
     - Renamed
   * - --cluster-status CLUSTER_STATUS
     -
                        Status command for cluster execution. This is only
                        considered in combination with the --cluster flag. If
                        provided, Snakemake will use the status command to
                        determine if a job has finished successfully or
                        failed. For this it is necessary that the submit
                        command provided to --cluster returns the cluster job
                        id. Then, the status command will be invoked with the
                        job id. Snakemake expects it to return 'success' if
                        the job was successful, 'failed' if the job failed and
                        'running' if the job still runs. (default: None)
     - --cluster-generic-status-cmd VALUE
       (Requires the cluster-generic plugin)
     -
                        Command for retrieving job status (default:
                        <dataclasses._MISSING_TYPE object at 0x7fc423088680>)
     - Renamed
   * - --cluster-cancel CLUSTER_CANCEL
     -
                        Specify a command that allows to stop currently
                        running jobs. The command will be passed a single
                        argument, the job id. (default: None)
     - --cluster-generic-cancel-cmd VALUE
       (Requires the cluster-generic plugin)
     -
                        Command for cancelling jobs. Expected to take one or
                        more jobids as arguments. (default:
                        <dataclasses._MISSING_TYPE object at 0x7fc423088680>)
     - Renamed
   * - --cluster-cancel-nargs CLUSTER_CANCEL_NARGS
     -
                        Specify maximal number of job ids to pass to
                        --cluster-cancel command, defaults to 1000. (default:
                        1000)
     - --cluster-generic-cancel-nargs VALUE
       (Requires the cluster-generic plugin)
     -
                        Number of jobids to pass to cancel_cmd. If more are
                        given, cancel_cmd will be called multiple times.
                        (default: <dataclasses._MISSING_TYPE object at
                        0x7fc423088680>)
     - Renamed
   * - --cluster-sidecar CLUSTER_SIDECAR
     -
                        Optional command to start a sidecar process during
                        cluster execution. Only active when --cluster is given
                        as well. (default: None)
     - --cluster-generic-sidecar-cmd VALUE
       (Requires the cluster-generic plugin)
     -
                        Command for sidecar process. (default:
                        <dataclasses._MISSING_TYPE object at 0x7fc423088680>)
     - Renamed


Profiles
^^^^^^^^

Profiles can now be versioned.
If your profile makes use of settings that are available in version 8 or later, use the filename ``config.v8+.yaml`` for the profile configuration (see :ref:`profiles <profiles>`).

API
^^^

The Snakemake API has been completely rewritten into a modern `dataclass <https://docs.python.org/3/library/dataclasses.html>`_ based approach.
The traditional central ``snakemake()`` function is gone.
For an example how to use the new API, check out the Snakemake CLI implementation `here <https://github.com/snakemake/snakemake/blob/04ec2c0262b2cb96cbcd7edbbb2596979c1703ae/snakemake/cli.py#L1767>`__.
