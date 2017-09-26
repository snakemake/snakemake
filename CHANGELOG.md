# Change Log

## [4.1.0] - 2017-09-26
### Added
- Support for configuration profiles. Profiles allow to specify default options, e.g., a cluster
  submission command. They can be used via 'snakemake --profile myprofile'. See the docs for details.
- GFAL remote provider. This allows to use GridFTP, SRM and any other protocol supported by GFAL for remote input and output files.
- Added --cluster-status flag that allows to specify a command that returns jobs status.
### Changed
- The scheduler now tries to get rid of the largest temp files first.
- The Docker image used for kubernetes support can now be configured at the command line.
- Rate-limiting for cluster interaction has be unified.
- S3 remote provider uses boto3.
- Resource functions can now use an additional `attempt` parameter, that contains the number of times this job has already been tried.
- Various minor fixes.

## [4.0.0] - 2017-07-24
### Added
- Cloud computing support via Kubernetes. Snakemake workflows can be executed transparently
  in the cloud, while storing input and output files within the cloud storage
  (e.g. S3 or Google Storage). I.e., this feature does not need a shared filesystem
  between the cloud notes, and thereby makes the setup really simple.
- WebDAV remote file support: Snakemake can now read and write from WebDAV. Hence,
  it can now, e.g., interact with Nextcloud or Owncloud.
- Support for default remote providers: define a remote provider to implicitly
  use for all input and output files.
- Added an option to only create conda environments instead of executing the workflow.
### Changed
- The number of files used for the metadata tracking of Snakemake (e.g., code, params, input changes) in the .snakemake directory has been reduced by a factor of 10, which should help with NFS and IO bottlenecks. This is a breaking change in the sense that Snakemake 4.x won't see the metadata of workflows executed with Snakemake 3.x. However, old metadata won't be overwritten, so that you can always go back and check things by installing an older version of Snakemake again.
- The google storage (GS) remote provider has been changed to use the google SDK.
  This is a breaking change, since the remote provider invocation has been simplified (see docs).
- Due to WebDAV support (which uses asyncio), Snakemake now requires Python 3.5 at least.
- Various minor bug fixes (e.g. for dynamic output files).


## [3.13.3] - 2017-06-23
### Changed
- Fix a followup bug in Namedlist where a single item was not returned as string.


## [3.13.2] - 2017-06-20
### Changed
- The --wrapper-prefix flag now also affects where the corresponding environment definition is fetched from.
- Fix bug where empty output file list was recognized as containing duplicates (issue #574).


## [3.13.1] - 2017-06-20
### Changed
- Fix --conda-prefix to be passed to all jobs.
- Fix cleanup issue with scripts that fail to download.

## [3.13.0] - 2017-06-12
### Added
- An NCBI remote provider. By this, you can seamlessly integrate any NCBI resouce (reference genome, gene/protein sequences, ...) as input file.
### Changed
- Snakemake now detects if automatically generated conda environments have to be recreated because the workflow has been moved to a new path.
- Remote functionality has been made more robust, in particular to avoid race conditions.
- `--config` parameter evaluation has been fixed for non-string types.
- The Snakemake docker container is now based on the official debian image.

## [3.12.0] - 2017-05-09
### Added
- Support for RMarkdown (.Rmd) in script directives.
- New option --debug-dag that prints all decisions while building the DAG of jobs. This helps to debug problems like cycles or unexpected MissingInputExceptions.
- New option --conda-prefix to specify the place where conda environments are stored.

### Changed
- Benchmark files now also include the maximal RSS and VMS size of the Snakemake process and all sub processes.
- Speedup conda environment creation.
- Allow specification, of DRMAA log dir.
- Pass cluster config to subworkflow.


## [3.11.2] - 2017-03-15
### Changed
- Fixed fix handling of local URIs with the wrapper directive.


## [3.11.1] - 2017-03-14
### Changed
- --touch ignores missing files
- Fixed handling of local URIs with the wrapper directive.


## [3.11.0] - 2017-03-08
### Added
- Param functions can now also refer to threads.
### Changed
- Improved tutorial and docs.
- Made conda integration more robust.
- None is converted to NULL in R scripts.


## [3.10.2] - 2017-02-28
### Changed
- Improved config file handling and merging.
- Output files can be referred in params functions (i.e. lambda wildcards, output: ...)
- Improved conda-environment creation.
- Jobs are cached, leading to reduced memory footprint.
- Fixed subworkflow handling in input functions.

## [3.10.0] - 2017-01-18
### Added
- Workflows can now be archived to a tarball with `snakemake --archive my-workflow.tar.gz`. The archive contains all input files, source code versioned with git and all software packages that are defined via conda environments. Hence, the archive allows to fully reproduce a workflow on a different machine. Such an archive can be uploaded to Zenodo, such that your workflow is secured in a self-contained, executable way for the future.
### Changed
- Improved logging.
- Reduced memory footprint.
- Added a flag to automatically unpack the output of input functions.
- Improved handling of HTTP redirects with remote files.
- Improved exception handling with DRMAA.
- Scripts referred by the script directive can now use locally defined external python modules.


## [3.9.1] - 2016-12-23
### Added
- Jobs can be restarted upon failure (--restart-times).
### Changed
- The docs have been restructured and improved. Now available under snakemake.readthedocs.org.
- Changes in scripts show up with --list-code-changes.
- Duplicate output files now cause an error.
- Various bug fixes.


## [3.9.0] - 2016-11-15
### Added
- Ability to define isolated conda software environments (YAML) per rule. Environment will be deployed by Snakemake upon workflow execution.
- Command line argument --wrapper-prefix in order to overwrite the default URL for looking up wrapper scripts.
### Changed
- --summary now displays the log files correspoding to each output file.
- Fixed hangups when using run directive and a large number of jobs
- Fixed pickling errors with anonymous rules and run directive.
- Various small bug fixes

## [3.8.2] - 2016-09-23
### Changed
- Add missing import in rules.py.
- Use threading only in cluster jobs.

## [3.8.1] - 2016-09-14
### Changed
- Snakemake now warns when using relative paths starting with "./".
- The option -R now also accepts an empty list of arguments.
- Bug fix when handling benchmark directive.
- Jobscripts exit with code 1 in case of failure. This should improve the error messages of cluster system.
- Fixed a bug in SFTP remote provider.


## [3.8.0] - 2016-08-26
### Added
- Wildcards can now be constrained by rule and globally via the new `wildcard_constraints` directive (see the [docs](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-wildcards)).
- Subworkflows now allow to overwrite their config file via the configfile directive in the calling Snakefile.
- A method `log_fmt_shell` in the snakemake proxy object that is available in scripts and wrappers allows to obtain a formatted string to redirect logging output from STDOUT or STDERR.
- Functions given to resources can now optionally contain an additional argument `input` that refers to the input files.
- Functions given to params can now optionally contain additional arguments `input` (see above) and `resources`. The latter refers to the resources.
- It is now possible to let items in shell commands be automatically quoted (see the [docs](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-rules)). This is usefull when dealing with filenames that contain whitespaces.
### Changed
- Snakemake now deletes output files before job exection. Further, it touches output files after job execution. This solves various problems with slow NFS filesystems.
- A bug was fixed that caused dynamic output rules to be executed multiple times when forcing their execution with -R.
- A bug causing double uploads with remote files was fixed. Various additional bug fixes related to remote files.
- Various minor bug fixes.

## [3.7.1] - 2016-05-16
### Changed
- Fixed a missing import of the multiprocessing module.

## [3.7.0] - 2016-05-05
### Added
- The entries in `resources` and the `threads` job attribute can now be callables that must return `int` values.
- Multiple `--cluster-config` arguments can be given to the Snakemake command line. Later one override earlier ones.
- In the API, multiple `cluster_config` paths can be given as a list, alternatively to the previous behaviour of expecting one string for this parameter.
- When submitting cluster jobs (either through `--cluster` or `--drmaa`), you can now use `--max-jobs-per-second` to limit the number of jobs being submitted (also available through Snakemake API). Some cluster installations have problems with too many jobs per second.
- Wildcard values are now printed upon job execution in addition to input and output files.
### Changed
- Fixed a bug with HTTP remote providers.

## [3.6.1] - 2016-04-08
### Changed
- Work around missing RecursionError in Python < 3.5
- Improved conversion of numpy and pandas data structures to R scripts.
- Fixed locking of working directory.

## [3.6.0] - 2016-03-10
### Added
- onstart handler, that allows to add code that shall be only executed before the actual workflow execution (not on dryrun).
- Parameters defined in the cluster config file are now accessible in the job properties under the key "cluster".
- The wrapper directive can be considered stable.
### Changed
- Allow to use rule/job parameters with braces notation in cluster config.
- Show a proper error message in case of recursion errors.
- Remove non-empty temp dirs.
- Don't set the process group of Snakemake in order to allow kill signals from parent processes to be propagated.
- Fixed various corner case bugs.
- The params directive no longer converts a list ``l`` implicitly to ``" ".join(l)``.

## [3.5.5] - 2016-01-23
### Added
- New experimental wrapper directive, which allows to refer to re-usable [wrapper scripts](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-wrappers). Wrappers are provided in the [Snakemake Wrapper Repository](https://bitbucket.org/snakemake/snakemake-wrappers).
- David Koppstein implemented two new command line options to constrain the execution of the DAG of job to sub-DAGs (--until and --omit-from).
### Changed
- Fixed various bugs, e.g. with shadow jobs and --latency-wait.

## [3.5.4] - 2015-12-04
### Changed
- The params directive now fully supports non-string parameters. Several bugs in the remote support were fixed.

## [3.5.3] - 2015-11-24
### Changed
- The missing remote module was added to the package.

## [3.5.2] - 2015-11-24
### Added
- Support for easy integration of external R and Python scripts via the new [script directive](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-external-scripts).
- Chris Tomkins-Tinch has implemented support for remote files: Snakemake can now handle input and output files from Amazon S3, Google Storage, FTP, SFTP, HTTP and Dropbox.
- Simon Ye has implemented support for sandboxing jobs with [shadow rules](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-shadow-rules).
### Changed
- Manuel Holtgrewe has fixed dynamic output files in combination with mutliple wildcards.
- It is now possible to add suffixes to all shell commands with shell.suffix("mysuffix").
- Job execution has been refactored to spawn processes only when necessary, resolving several problems in combination with huge workflows consisting of thousands of jobs and reducing the memory footprint.
- In order to reflect the new collaborative development model, Snakemake has moved from my personal bitbucket account to http://snakemake.bitbucket.org.

## [3.4.2] - 2015-09-12
### Changed
- Willem Ligtenberg has reduced the memory usage of Snakemake.
- Per Unneberg has improved config file handling to provide a more intuitive overwrite behavior.
- Simon Ye has improved the test suite of Snakemake and helped with setting up continuous integration via Codeship.
- The cluster implementation has been rewritten to use only a single thread to wait for jobs. This avoids failures with large numbers of jobs.
- Benchmarks are now writing tab-delimited text files instead of JSON.
- Snakemake now always requires to set the number of jobs with -j when in cluster mode. Set this to a high value if your cluster does not have restrictions.
- The Snakemake Conda package has been moved to the bioconda channel.
- The handling of Symlinks was improved, which made a switch to Python 3.3 as the minimum required Python version necessary.

## [3.4.1] - 2015-08-05
### Changed
- This release fixes a bug that caused named input or output files to always be returned as lists instead of single files.

## [3.4] - 2015-07-18
### Added
- This release adds support for executing jobs on clusters in synchronous mode (e.g. qsub -sync). Thanks to David Alexander for implementing this.
- There is now vim syntax highlighting support (thanks to Jay Hesselberth).
- Snakemake is now available as Conda package.
### Changed
- Lots of bugs have been fixed. Thanks go to e.g. David Koppstein, Marcel Martin, John Huddleston and Tao Wen for helping with useful reports and debugging.

See [here](https://bitbucket.org/snakemake/snakemake/wiki/News-Archive) for older changes.
