__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from collections import defaultdict
import os
import re
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path
from typing import List, Mapping, Optional, Set, Union
from snakemake import caching
from snakemake_interface_executor_plugins.settings import ExecMode
from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry
from snakemake_interface_executor_plugins.utils import is_quoted, maybe_base64
from snakemake_interface_storage_plugins.registry import StoragePluginRegistry
from snakemake_interface_report_plugins.registry import ReportPluginRegistry

from snakemake_interface_logger_plugins.registry import LoggerPluginRegistry


import snakemake.common.argparse
from snakemake import logging
from snakemake.api import (
    SnakemakeApi,
    resolve_snakefile,
)
from snakemake.common import (
    SNAKEFILE_CHOICES,
    __version__,
    async_run,
    get_appdirs,
    get_container_image,
    parse_key_value_arg,
)
from snakemake.exceptions import (
    CliException,
    ResourceScopesException,
    print_exception,
)
from snakemake.resources import (
    DefaultResources,
    ParsedResource,
    ResourceScopes,
    eval_resource_expression,
    parse_resources,
)
from snakemake.settings.types import (
    Batch,
    ChangeType,
    ConfigSettings,
    DAGSettings,
    DeploymentMethod,
    DeploymentSettings,
    ExecutionSettings,
    GroupSettings,
    MaxJobsPerTimespan,
    NotebookEditMode,
    OutputSettings,
    PreemptibleRules,
    Quietness,
    RemoteExecutionSettings,
    RerunTrigger,
    ResourceSettings,
    SchedulingSettings,
    SharedFSUsage,
    StorageSettings,
    WorkflowSettings,
    StrictDagEvaluation,
    PrintDag,
)
from snakemake.target_jobs import parse_target_jobs_cli_args
from snakemake.utils import available_cpu_count, update_config


def parse_size_in_bytes(value):
    from humanfriendly import parse_size

    return parse_size(value)


def parse_timespan(value):
    from humanfriendly import parse_timespan

    return parse_timespan(value)


def expandvars(atype):
    def inner(args):
        if isinstance(args, list):
            return atype([os.path.expandvars(arg) for arg in args])
        elif isinstance(args, str):
            return atype(os.path.expandvars(args))
        elif args is None:
            return None
        else:
            return atype(args)

    return inner


def optional_str(arg):
    if arg is None or arg == "none":
        return None
    else:
        return arg


def parse_set_threads(args):
    def fallback(orig_value):
        value = eval_resource_expression(orig_value, threads_arg=False)
        return ParsedResource(value=value, orig_arg=orig_value)

    return parse_set_ints(
        args,
        "Invalid threads definition: entries have to be defined as RULE=THREADS pairs "
        "(with THREADS being a positive integer).",
        fallback=fallback,
    )


def parse_consider_ancient(
    args: Optional[List[str]],
) -> Mapping[str, Set[Union[str, int]]]:
    """Parse command line arguments for marking input files as ancient.

    Args:
        args: List of RULE=INPUTITEMS pairs, where INPUTITEMS is a comma-separated list
              of input item names or indices (0-based).

    Returns:
        A mapping of rules to sets of their ancient input items.

    Raises:
        ValueError: If the format is invalid or values cannot be parsed.
    """
    errmsg = (
        "Invalid --consider-ancient definition: entries have to be defined as "
        "RULE=INPUTITEMS pairs, with INPUTITEMS being a list of input items of the "
        "rule (given as name or index (0-based)), separated by commas."
    )

    def parse_item(item: str) -> Union[str, int]:
        try:
            return int(item)
        except ValueError:
            if item.isidentifier():
                return item
            else:
                raise ValueError(f"{errmsg} (Unparsable value: {repr(item)})")

    consider_ancient = defaultdict(set)

    if args is not None:
        for entry in args:
            rule, items = parse_key_value_arg(entry, errmsg=errmsg, strip_quotes=True)
            if not rule.isidentifier():
                raise ValueError(f"{errmsg} (Invalid rule name: {repr(rule)})")
            items = items.split(",")
            consider_ancient[rule] = {parse_item(item) for item in items}
    return consider_ancient


def parse_set_resources(args):
    errmsg = (
        "Invalid resource definition: entries have to be defined as RULE:RESOURCE=VALUE, with "
        "VALUE being a positive integer a quoted string, or a Python expression (e.g. min(max(2*input.size_mb, 1000), 8000))."
    )

    from collections import defaultdict

    assignments = defaultdict(dict)
    if args is not None:
        for entry in args:
            key, orig_value = parse_key_value_arg(
                entry, errmsg=errmsg, strip_quotes=False
            )
            key = key.split(":")
            if len(key) != 2:
                raise ValueError(errmsg)
            rule, resource = key
            if is_quoted(orig_value):
                # value is a string, just keep it but remove surrounding quotes
                value = orig_value[1:-1]
            else:
                try:
                    value = int(orig_value)
                except ValueError:
                    value = eval_resource_expression(orig_value)
            if isinstance(value, int) and value < 0:
                raise ValueError(errmsg)
            assignments[rule][resource] = ParsedResource(
                value=value, orig_arg=orig_value
            )
    return assignments


def parse_set_scatter(args):
    return parse_set_ints(
        args,
        "Invalid scatter definition: entries have to be defined as NAME=SCATTERITEMS pairs "
        "(with SCATTERITEMS being a positive integer).",
    )


def parse_set_resource_scope(args):
    err_msg = (
        "Invalid resource scopes: entries must be defined as RESOURCE=SCOPE pairs, "
        "where SCOPE is either 'local', 'global', or 'excluded'"
    )
    if args is not None:
        try:
            return ResourceScopes(
                parse_key_value_arg(entry, errmsg=err_msg) for entry in args
            )
        except ResourceScopesException as err:
            invalid_resources = ", ".join(
                f"'{res}={scope}'" for res, scope in err.invalid_resources.items()
            )
            raise ValueError(f"{err.msg} (got {invalid_resources})")

    return ResourceScopes()


def parse_set_ints(arg, errmsg, fallback=None):
    assignments = dict()
    if arg is not None:
        for entry in arg:
            key, value = parse_key_value_arg(entry, errmsg=errmsg)
            try:
                value = int(value)
            except ValueError:
                if fallback is not None:
                    try:
                        value = fallback(value)
                    except Exception as e:
                        raise ValueError(f"{errmsg} Cause: {e}")
                else:
                    raise ValueError(errmsg)
            if isinstance(value, int) and value < 0:
                raise ValueError(errmsg)
            assignments[key] = value
    return assignments


def parse_batch(arg):
    errmsg = "Invalid batch definition: batch entry has to be defined as RULE=BATCH/BATCHES (with integers BATCH <= BATCHES, BATCH >= 1)."
    if arg is not None:
        rule, batchdef = parse_key_value_arg(arg, errmsg=errmsg)
        try:
            batch, batches = batchdef.split("/")
            batch = int(batch)
            batches = int(batches)
        except ValueError:
            raise ValueError(errmsg)
        if batch > batches or batch < 1:
            raise ValueError(errmsg)
        return Batch(rule, batch, batches)
    return None


def parse_groups(args):
    errmsg = "Invalid groups definition: entries have to be defined as RULE=GROUP pairs"
    overwrite_groups = dict()
    if args is not None:
        for entry in args:
            rule, group = parse_key_value_arg(entry, errmsg=errmsg)
            overwrite_groups[rule] = group
    return overwrite_groups


def parse_group_components(args):
    errmsg = "Invalid group components definition: entries have to be defined as GROUP=COMPONENTS pairs (with COMPONENTS being a positive integer)"
    group_components = dict()
    if args is not None:
        for entry in args:
            group, count = parse_key_value_arg(entry, errmsg=errmsg)
            try:
                count = int(count)
            except ValueError:
                raise ValueError(errmsg)
            if count <= 0:
                raise ValueError(errmsg)
            group_components[group] = count
    return group_components


def _bool_parser(value):
    if value == "True":
        return True
    elif value == "False":
        return False
    raise ValueError


def parse_config(entries):
    """Parse config from args."""
    import yaml

    yaml_base_load = lambda s: yaml.load(s, Loader=yaml.loader.BaseLoader)
    parsers = [int, float, _bool_parser, yaml_base_load, str]
    config = dict()
    if entries:
        valid = re.compile(r"[a-zA-Z_][\w-]*\w$")
        for entry in entries:
            key, val = parse_key_value_arg(
                entry,
                errmsg="Invalid config definition: Config entries have to be defined as name=value pairs.",
            )
            if not valid.match(key):
                raise ValueError(
                    "Invalid config definition: Config entry must start with a valid identifier."
                )
            v = None
            if val == "":
                update_config(config, {key: v})
                continue
            for parser in parsers:
                try:
                    v = parser(val)
                    # avoid accidental interpretation as function
                    if not callable(v):
                        break
                except:
                    pass
            assert v is not None
            update_config(config, {key: v})
    return config


def parse_cores(cores):
    if cores == "all":
        return available_cpu_count()
    try:
        return int(cores)
    except ValueError:
        raise CliException(
            "Error parsing number of cores (--cores, -c): must be integer or 'all'."
        )


def parse_jobs(jobs):
    if jobs == "unlimited":
        return sys.maxsize
    try:
        return int(jobs)
    except ValueError:
        raise CliException(
            "Error parsing number of jobs (--jobs, -j): must be integer."
        )


def get_profile_dir(profile: str) -> (Path, Path):
    config_pattern = re.compile(r"config(.v(?P<min_major>\d+)\+)?.yaml")

    def get_config_min_major(filename):
        m = config_pattern.match(filename)
        if m:
            min_major = m.group("min_major")
            if min_major is None:
                return 0
            min_major = int(min_major)

            return min_major
        return None

    dirs = get_appdirs()
    if os.path.exists(profile):
        parent_dir = os.path.dirname(profile) or "."
        search_dirs = [parent_dir]
        profile = os.path.basename(profile)
    else:
        search_dirs = [os.getcwd(), dirs.user_config_dir, dirs.site_config_dir]
    for d in search_dirs:
        profile_candidate = Path(d) / profile
        if profile_candidate.exists():
            files = os.listdir(profile_candidate)
            # If versioneer cannot get the real version it will return something
            # like "0+untagged.5410.g40ffe59" - this should only occur in testing scenarios
            curr_major = int(__version__.split(".")[0].split("+")[0])
            config_files = {
                f: min_major
                for f, min_major in zip(files, map(get_config_min_major, files))
                if min_major is not None and curr_major >= min_major
            }
            if config_files:
                config_file = max(config_files, key=config_files.get)
                return profile_candidate, profile_candidate / config_file


def get_argument_parser(profiles=None):
    """Generate and return argument parser."""
    from snakemake.profiles import ProfileConfigFileParser

    dirs = get_appdirs()
    config_files = []
    if profiles:
        for profile in profiles:
            if profile == "":
                print("Error: invalid profile name.", file=sys.stderr)
                exit(1)

            profile_entry = get_profile_dir(profile)
            if profile_entry is not None:
                _profile_dir, config_file = profile_entry
                config_files.append(config_file)
            else:
                print(
                    "Error: profile given but no config.yaml found. "
                    "Profile has to be given as either absolute path, relative "
                    "path or name of a directory available in either "
                    "{site} or {user}.".format(
                        site=dirs.site_config_dir, user=dirs.user_config_dir
                    ),
                    file=sys.stderr,
                )
                exit(1)

    parser = snakemake.common.argparse.ArgumentParser(
        description="Snakemake is a Python based language and execution "
        "environment for GNU Make-like workflows.",
        formatter_class=snakemake.common.argparse.ArgumentDefaultsHelpFormatter,
        default_config_files=config_files,
        config_file_parser_class=ProfileConfigFileParser,
    )

    group_exec = parser.add_argument_group("EXECUTION")

    group_exec.add_argument(
        "targets",
        nargs="*",
        default=set(),
        help="Targets to build. May be rules or files.",
    )

    group_exec.add_argument(
        "--dry-run",
        "--dryrun",
        "-n",
        dest="dryrun",
        action="store_true",
        help="Do not execute anything, and display what would be done. "
        "If you have a very large workflow, use `--dry-run --quiet` to just "
        "print a summary of the DAG of jobs.",
    )

    group_exec.add_argument(
        "--profile",
        help=f"Name of profile to use for configuring Snakemake. Snakemake will search for a corresponding folder in `{dirs.site_config_dir}` and `{dirs.user_config_dir}`. Alternatively, this can be an absolute or relative path. The profile folder has to contain a file `config.yaml`. This file can be used to set default values for command line options in YAML format. For example, `--cluster qsub` becomes `cluster: qsub` in the YAML file. Profiles can be obtained from https://github.com/snakemake-profiles. The profile can also be set via the environment variable `$SNAKEMAKE_PROFILE`. To override this variable and use no profile at all, provide the value `none` to this argument.",
        env_var="SNAKEMAKE_PROFILE",
    )

    group_exec.add_argument(
        "--workflow-profile",
        help="Path (relative to current directory) to workflow specific profile folder to use for configuring Snakemake with parameters specific for this workflow (like resources). If this flag is not used, Snakemake will by default use `profiles/default` if present (searched both relative to current directory and relative to Snakefile, in this order). For skipping any workflow specific profile provide the special value `none`. Settings made in the workflow profile will override settings made in the general profile (see `--profile`). The profile folder has to contain a file `config.yaml`. This file can be used to set default values for command line options in YAML format. For example, `--executor slurm` becomes `executor: slurm` in the YAML file. It is advisable to use the workflow profile to set or overwrite e.g. workflow specific resources like the amount of threads of a particular rule or the amount of memory needed. Note that in such cases, the arguments may be given as nested YAML mappings in the profile, e.g. `set-threads: myrule: 4` instead of `set-threads: myrule=4`.",
    )

    group_exec.add_argument(
        "--cache",
        nargs="*",
        metavar="RULE",
        help="Store output files of given rules in a central cache given by the environment "
        f"variable `${caching.LOCATION_ENVVAR}`. Likewise, retrieve output files of the given rules "
        "from this cache if they have been created before (by anybody writing to the same cache), "
        "instead of actually executing the rules. Output files are identified by hashing all "
        "steps, parameters and software stack (conda envs or containers) needed to create them.",
    )

    group_exec.add_argument(
        "--snakefile",
        "-s",
        metavar="FILE",
        type=Path,
        help=(
            "The workflow definition in form of a snakefile. "
            "Usually, you should not need to specify this. "
            "By default, Snakemake will search for {} "
            "beneath the current working "
            "directory, in this order. "
            "Only if you definitely want a different layout, "
            "you need to use this parameter."
        ).format(", ".join(map("`{}`".format, SNAKEFILE_CHOICES))),
    )
    group_exec.add_argument(
        "--cores",
        "-c",
        action="store",
        metavar="N",
        type=parse_cores,
        help=(
            "Use at most N CPU cores/jobs in parallel. "
            "If N is omitted or `all`, the limit is set to the number of "
            "available CPU cores. "
            "In case of cluster/cloud execution, this argument sets the maximum number "
            "of cores requested from the cluster or cloud scheduler. (See "
            "https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#"
            "resources-remote-execution for more info.) "
            "This number is available to rules via workflow.cores."
        ),
    )
    group_exec.add_argument(
        "--jobs",
        "-j",
        metavar="N",
        action="store",
        type=parse_jobs,
        help=(
            "Use at most N CPU cluster/cloud jobs in parallel. For local execution this is "
            "an alias for `--cores` (it is though recommended to use `--cores` in that case). "
            "Note: Set to `unlimited` to allow any number of parallel jobs."
        ),
    )
    group_exec.add_argument(
        "--local-cores",
        action="store",
        metavar="N",
        type=int,
        help=(
            "In cluster/cloud mode, use at most N cores of the host machine in parallel "
            "(default: number of CPU cores of the host). The cores are used to execute "
            "local rules. This option is ignored when not in cluster/cloud mode."
        ),
    )
    group_exec.add_argument(
        "--resources",
        "--res",
        nargs="+",
        metavar="NAME=INT",
        default=dict(),
        parse_func=parse_resources,
        help=(
            "Define additional resources that shall constrain the scheduling "
            "analogously to `--cores` (see above). A resource is defined as "
            "a name and an integer value. E.g. `--resources mem_mb=1000`. Rules can "
            "use resources by defining the resource keyword, e.g. "
            "`resources: mem_mb=600`. If now two rules require 600 of the resource "
            "`mem_mb` they won't be run in parallel by the scheduler. In "
            "cluster/cloud mode, this argument will also constrain the amount of "
            "resources requested from the server. (See "
            "https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#"
            "resources-remote-execution for more info.)"
        ),
    )
    group_exec.add_argument(
        "--set-threads",
        metavar="RULE=THREADS",
        nargs="+",
        default=dict(),
        parse_func=maybe_base64(parse_set_threads),
        help="Overwrite thread usage of rules. This allows to fine-tune workflow "
        "parallelization. In particular, this is helpful to target certain cluster nodes "
        "by e.g. shifting a rule to use more, or less threads than defined in the workflow. "
        "Thereby, THREADS has to be a positive integer, and RULE has to be the name of the rule.",
    )
    group_exec.add_argument(
        "--max-threads",
        type=int,
        help="Define a global maximum number of threads available to any rule. Rules "
        "requesting more threads (via the threads keyword) will have their values "
        "reduced to the maximum. This can be useful when you want to restrict the "
        "maximum number of threads without modifying the workflow definition or "
        "overwriting rules individually with `--set-threads`.",
    )
    group_exec.add_argument(
        "--set-resources",
        metavar="RULE:RESOURCE=VALUE",
        nargs="+",
        default=dict(),
        parse_func=maybe_base64(parse_set_resources),
        help="Overwrite resource usage of rules. This allows to fine-tune workflow "
        "resources. In particular, this is helpful to target certain cluster nodes "
        "by e.g. defining a certain partition for a rule, or overriding a temporary directory. "
        "Thereby, VALUE has to be a positive integer or a string, RULE has to be the name of the "
        "rule, and RESOURCE has to be the name of the resource.",
    )
    group_exec.add_argument(
        "--set-scatter",
        metavar="NAME=SCATTERITEMS",
        nargs="+",
        default=dict(),
        parse_func=parse_set_scatter,
        help="Overwrite number of scatter items of scattergather processes. This allows to fine-tune "
        "workflow parallelization. Thereby, SCATTERITEMS has to be a positive integer, and NAME has to be "
        "the name of the scattergather process defined via a scattergather directive in the workflow.",
    )
    group_exec.add_argument(
        "--set-resource-scopes",
        metavar="RESOURCE=[global|local]",
        nargs="+",
        default=dict(),
        parse_func=parse_set_resource_scope,
        help="Overwrite resource scopes. A scope determines how a constraint is "
        "reckoned in cluster execution. With RESOURCE=local, a constraint applied to "
        "RESOURCE using `--resources` will be considered the limit for each group "
        "submission. With RESOURCE=global, the constraint will apply across all groups "
        "cumulatively. By default, only `mem_mb` and `disk_mb` are considered local, "
        "all other resources are global. This may be modified in the snakefile using "
        "the `resource_scopes:` directive. Note that number of threads, specified via "
        "`--cores`, is always considered local. (See "
        "https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#"
        "resources-remote-execution for more info)",
    )
    group_exec.add_argument(
        "--default-resources",
        "--default-res",
        nargs="*",
        metavar="NAME=INT",
        parse_func=maybe_base64(DefaultResources),
        help=(
            "Define default values of resources for rules that do not define their own values. "
            "In addition to plain integers, python expressions over inputsize are allowed (e.g. `2*input.size_mb`). "
            "The inputsize is the sum of the sizes of all input files of a rule. "
            "By default, Snakemake assumes a default for mem_mb, disk_mb, and tmpdir (see below). "
            "This option allows to add further defaults (e.g. account and partition for slurm) or to overwrite these default values. "
            "The defaults are `mem_mb=min(max(2*input.size_mb, 1000), 8000)`, `disk_mb=max(2*input.size_mb, 1000)` "
            "(i.e., default disk and mem usage is twice the input file size but at least 1GB), and "
            "the system temporary directory (as given by $TMPDIR, $TEMP, or $TMP) is used for the tmpdir resource. "
            "The tmpdir resource is automatically used by shell commands, scripts and wrappers to store temporary data (as it is "
            "mirrored into $TMPDIR, $TEMP, and $TMP for the executed subprocesses). "
            "If this argument is not specified at all, Snakemake just uses the tmpdir resource as outlined above. "
            "The tmpdir resource can also be overwritten in the same way as e.g. mem_mb above. "
            "Thereby, it is even possible to use shutil.disk_usage(system_tmpdir).free and comparing this to input.size in order to "
            "determine if one can expect the system_tmpdir to be big enough and switch to another tmpdir in case it is not. "
        ),
    )

    group_exec.add_argument(
        "--preemptible-rules",
        nargs="*",
        parse_func=set,
        help=(
            "Define which rules shall use a preemptible machine which can be prematurely killed by e.g. a cloud provider (also called spot instances). "
            "This is currently only supported by the Google Life Sciences executor and ignored by all other executors. "
            "If no rule names are provided, all rules are considered to be preemptible. "
        ),
    )

    group_exec.add_argument(
        "--preemptible-retries",
        type=int,
        help="Number of retries that shall be made in order to finish a job from of rule that has been marked as preemptible via the --preemptible-rules setting.",
    )
    group_exec.add_argument(
        "--configfile",
        "--configfiles",
        nargs="+",
        metavar="FILE",
        default=list(),
        type=Path,
        help=(
            "Specify or overwrite the config file of the workflow (see the docs). "
            "Values specified in JSON or YAML format are available in the global config "
            "dictionary inside the workflow. Multiple files overwrite each other in "
            "the given order. Thereby missing keys in previous config files are extended by "
            "following configfiles. Note that this order also includes a config file defined "
            "in the workflow definition itself (which will come first)."
        ),
    )
    group_exec.add_argument(
        "--config",
        "-C",
        nargs="*",
        metavar="KEY=VALUE",
        help=(
            "Set or overwrite values in the workflow config object. "
            "The workflow config object is accessible as variable config inside "
            "the workflow. Default values can be set by providing a YAML JSON file "
            "(see `--configfile` and Documentation)."
        ),
    )
    group_exec.add_argument(
        "--replace-workflow-config",
        action="store_true",
        help=(
            "Config files provided via command line do not update and extend the config "
            "dictionary of the workflow but instead fully replace it. Keys that are not "
            "defined in the provided config files will be undefined even if specified "
            "within the workflow config."
        ),
    )
    group_exec.add_argument(
        "--envvars",
        nargs="+",
        metavar="VARNAME",
        parse_func=set,
        default=set(),
        help="Environment variables to pass to cloud jobs.",
    )
    group_exec.add_argument(
        "--directory",
        "-d",
        metavar="DIR",
        type=Path,
        help=(
            "Specify working directory (relative paths in "
            "the snakefile will use this as their origin)."
        ),
    )
    group_exec.add_argument(
        "--touch",
        "-t",
        action="store_true",
        help=(
            "Touch output files (mark them up to date without really "
            "changing them) instead of running their commands. This is "
            "used to pretend that the rules were executed, in order to "
            "fool future invocations of snakemake. Fails if a file does "
            "not yet exist. Note that this will only touch files that would "
            "otherwise be recreated by Snakemake (e.g. because their input "
            "files are newer). For enforcing a touch, combine this with "
            "`--force`, `--forceall`, or `--forcerun`. Note however that you lose "
            "the provenance information when the files have been created in "
            "reality. Hence, this should be used only as a last resort."
        ),
    )
    group_exec.add_argument(
        "--keep-going",
        "-k",
        action="store_true",
        help="Go on with independent jobs if a job fails.",
    )
    group_exec.add_argument(
        "--rerun-triggers",
        nargs="+",
        choices=RerunTrigger.choices(),
        default=RerunTrigger.all(),
        parse_func=RerunTrigger.parse_choices_set,
        help="Define what triggers the rerunning of a job. By default, "
        "all triggers are used, which guarantees that results are "
        "consistent with the workflow code and configuration. If you "
        "rather prefer the traditional way of just considering "
        "file modification dates, use `--rerun-trigger mtime`.",
    )
    group_exec.add_argument(
        "--force",
        "-f",
        action="store_true",
        help=(
            "Force the execution of the selected target or the first rule "
            "regardless of already created output."
        ),
    )
    group_exec.add_argument(
        "--executor",
        "-e",
        help="Specify a custom executor, available via an executor plugin: snakemake_executor_<name>",
        choices=ExecutorPluginRegistry().plugins.keys(),
    )
    group_exec.add_argument(
        "--forceall",
        "-F",
        action="store_true",
        help=(
            "Force the execution of the selected (or the first) rule and "
            "all rules it is dependent on regardless of already created "
            "output."
        ),
    )
    group_exec.add_argument(
        "--forcerun",
        "-R",
        nargs="*",
        metavar="TARGET",
        parse_func=set,
        default=set(),
        help=(
            "Force the re-execution or creation of the given rules or files."
            " Use this option if you changed a rule and want to have all its "
            "output in your workflow updated."
        ),
    )
    group_exec.add_argument(
        "--consider-ancient",
        metavar="RULE=INPUTITEMS",
        nargs="+",
        default=dict(),
        parse_func=parse_consider_ancient,
        help="Consider given input items of given rules as ancient, i.e. not triggering "
        "re-runs if they are newer than the output files. "
        "Putting this into a workflow specific profile (or specifying as argument) "
        "allows to overrule rerun triggers caused by file modification dates where the "
        "user knows better. RULE is the name of the rule, INPUTITEMS is a comma "
        "separated list of input items of the rule (given as name or index (0-based)).",
    )

    group_exec.add_argument(
        "--prioritize",
        "-P",
        nargs="+",
        metavar="TARGET",
        parse_func=set,
        default=set(),
        help=(
            "Tell the scheduler to assign creation of given targets "
            "(and all their dependencies) highest priority."
        ),
    )
    group_exec.add_argument(
        "--batch",
        metavar="RULE=BATCH/BATCHES",
        type=parse_batch,
        help=(
            "Only create the given BATCH of the input files of the given RULE. "
            "This can be used to iteratively run parts of very large workflows. "
            "Only the execution plan of the relevant part of the workflow has to "
            "be calculated, thereby speeding up DAG computation. "
            "It is recommended to provide the most suitable rule for batching when "
            "documenting a workflow. It should be some aggregating rule that "
            "would be executed only once, and has a large number of input files. "
            "For example, it can be a rule that aggregates over samples."
        ),
    )
    group_exec.add_argument(
        "--until",
        "-U",
        nargs="+",
        metavar="TARGET",
        parse_func=set,
        default=set(),
        help=(
            "Runs the pipeline until it reaches the specified rules or "
            "files. Only runs jobs that are dependencies of the specified "
            "rule or files, does not run sibling DAGs. "
        ),
    )
    group_exec.add_argument(
        "--omit-from",
        "-O",
        nargs="+",
        metavar="TARGET",
        parse_func=set,
        default=set(),
        help=(
            "Prevent the execution or creation of the given rules or files "
            "as well as any rules or files that are downstream of these targets "
            "in the DAG. Also runs jobs in sibling DAGs that are independent of the "
            "rules or files specified here."
        ),
    )
    group_exec.add_argument(
        "--rerun-incomplete",
        "--ri",
        action="store_true",
        help=("Re-run all jobs the output of which is recognized as incomplete."),
    )
    group_exec.add_argument(
        "--shadow-prefix",
        metavar="DIR",
        help=(
            "Specify a directory in which the `shadow` directory is created. "
            "If not supplied, the value is set to the `.snakemake` directory relative "
            "to the working directory."
        ),
    )

    group_exec.add_argument(
        "--strict-dag-evaluation",
        nargs="+",
        choices=StrictDagEvaluation.choices(),
        default=set(),
        parse_func=StrictDagEvaluation.parse_choices_set,
        help="Strict evaluation of rules' correctness even when not required to produce the output files. ",
    )

    try:
        import pulp

        lp_solvers = pulp.listSolvers(onlyAvailable=True)
    except ImportError:
        # Dummy list for the case that pulp is not available
        # This only happened when building docs.
        lp_solvers = ["COIN_CMD"]
    recommended_lp_solver = "COIN_CMD"

    group_exec.add_argument(
        "--scheduler",
        default="greedy" if recommended_lp_solver not in lp_solvers else "ilp",
        nargs="?",
        choices=["ilp", "greedy"],
        help=(
            "Specifies if jobs are selected by a greedy algorithm or by solving an ilp. "
            "The ilp scheduler aims to reduce runtime and hdd usage by best possible use of resources."
        ),
    )

    group_exec.add_argument(
        "--scheduler-ilp-solver",
        default=recommended_lp_solver,
        choices=lp_solvers,
        help=("Specifies solver to be utilized when selecting ilp-scheduler."),
    )

    group_exec.add_argument(
        "--conda-base-path",
        help="Path of conda base installation (home of conda, mamba, activate) (internal use only).",
    )

    group_exec.add_argument(
        "--no-subworkflows",
        "--nosw",
        action="store_true",
        help=("Do not evaluate or execute subworkflows."),
    )
    group_exec.add_argument(
        "--precommand",
        help="Only used in case of remote execution. Command to be executed before "
        "Snakemake executes each job on the remote compute node.",
    )

    group_group = parser.add_argument_group("GROUPING")
    group_group.add_argument(
        "--groups",
        nargs="+",
        parse_func=parse_groups,
        default=dict(),
        help="Assign rules to groups (this overwrites any "
        "group definitions from the workflow).",
    )
    group_group.add_argument(
        "--group-components",
        nargs="+",
        parse_func=parse_group_components,
        default=dict(),
        help="Set the number of connected components a group is "
        "allowed to span. By default, this is 1, but this flag "
        "allows to extend this. This can be used to run e.g. 3 "
        "jobs of the same rule in the same group, although they "
        "are not connected. It can be helpful for putting together "
        "many small jobs or benefitting of shared memory setups.",
    )
    group_report = parser.add_argument_group("REPORTS")

    group_report.add_argument(
        "--report",
        nargs="?",
        const="report.html",
        metavar="FILE",
        type=Path,
        help="Create a self-contained HTML report with default statistics, "
        "provenance information and user-specified results. "
        "For smaller datasets with a limited report complexity, you can specify "
        "an `.html` file and all results will be embedded directly into this file. "
        "For customized reports on larger sample sizes, it makes more sense to "
        "specify a `.zip` file. The resulting archive will spread the contents "
        "across a folder structure, for a quicker loading of individual results. "
        "You can unpack this archive anywhere and open the `report.html` file in "
        "its main folder to view the report in any web browser.",
    )
    group_report.add_argument(
        "--report-after-run",
        action="store_true",
        help="After finishing the workflow, directly create the report. "
        "It is required to provide --report.",
    )
    group_report.add_argument(
        "--report-stylesheet",
        metavar="CSSFILE",
        type=Path,
        help="Custom stylesheet to use for report. In particular, this can be used for "
        "branding the report with e.g. a custom logo, see docs.",
    )
    group_report.add_argument(
        "--reporter",
        metavar="PLUGIN",
        help="Specify a custom report plugin. By default, Snakemake's builtin html "
        "reporter will be used. For custom reporters, check out their command line "
        "options starting with `--report-`.",
    )

    group_notebooks = parser.add_argument_group("NOTEBOOKS")

    group_notebooks.add_argument(
        "--draft-notebook",
        metavar="TARGET",
        help="Draft a skeleton notebook for the rule used to generate the given target file. This notebook "
        "can then be opened in a jupyter server, executed and implemented until ready. After saving, it "
        "will automatically be reused in non-interactive mode by Snakemake for subsequent jobs.",
    )
    group_notebooks.add_argument(
        "--edit-notebook",
        metavar="TARGET",
        help="Interactively edit the notebook associated with the rule used to generate the given target file. "
        "This will start a local jupyter notebook server. "
        "Any changes to the notebook should be saved, and the server has to be stopped by "
        "closing the notebook and hitting the `Quit` button on the jupyter dashboard. "
        "Afterwards, the updated notebook will be automatically stored in the path defined in the rule. "
        "If the notebook is not yet present, this will create an empty draft. ",
    )
    group_notebooks.add_argument(
        "--notebook-listen",
        metavar="IP:PORT",
        default="localhost:8888",
        help="The IP address and PORT the notebook server used for editing the notebook (`--edit-notebook`) will listen on.",
    )

    group_utils = parser.add_argument_group("UTILITIES")
    group_utils.add_argument(
        "--lint",
        nargs="?",
        const="text",
        choices=["text", "json"],
        help="Perform linting on the given workflow. This will print snakemake "
        "specific suggestions to improve code quality (work in progress, more lints "
        "to be added in the future). If no argument is provided, plain text output is used.",
    )
    group_utils.add_argument(
        "--generate-unit-tests",
        nargs="?",
        const=".tests/unit",
        metavar="TESTPATH",
        type=Path,
        help="Automatically generate unit tests for each workflow rule. "
        "This assumes that all input files of each job are already present. "
        "Jobs without present input files will be skipped (a warning will be issued). "
        "For each rule, one test case will be created and, after "
        "successful execution, tests can be run with `pytest TESTPATH`.",
    )
    group_utils.add_argument(
        "--containerize",
        action="store_true",
        help="Print a Dockerfile that provides an execution environment for the workflow, including all "
        "conda environments.",
    )
    group_utils.add_argument(
        "--export-cwl",
        action="store",
        metavar="FILE",
        help="Compile workflow to CWL and store it in given FILE.",
    )
    group_utils.add_argument(
        "--list-rules",
        "--list",
        "-l",
        action="store_true",
        help="Show available rules in given Snakefile.",
    )
    group_utils.add_argument(
        "--list-target-rules",
        "--lt",
        action="store_true",
        help="Show available target rules in given Snakefile.",
    )
    group_exec.add_argument(
        "--dag",
        nargs="?",
        choices=PrintDag.choices(),
        const=str(PrintDag.DOT),
        default=None,
        help="Do not execute anything and print the directed "
        "acyclic graph of jobs in the dot language or in mermaid-js. Recommended "
        "use on Unix systems: `snakemake --dag | dot | display`. "
        "Note print statements in your Snakefile may interfere "
        "with visualization.",
    )
    group_utils.add_argument(
        "--rulegraph",
        nargs="?",
        choices=PrintDag.choices(),
        const=str(PrintDag.DOT),
        default=None,
        help="Do not execute anything and print the dependency graph "
        "of rules in the dot language or in mermaid-js. This will be less "
        "crowded than above DAG of jobs, but also show less information. "
        "Note that each rule is displayed once, hence the displayed graph will be "
        "cyclic if a rule appears in several steps of the workflow. "
        "Use this if above option leads to a DAG that is too large. "
        "Recommended use on Unix systems: snakemake `--rulegraph | dot | display`. "
        "Note print statements in your Snakefile may interfere "
        "with visualization.",
    )
    group_utils.add_argument(
        "--filegraph",
        action="store_true",
        help="Do not execute anything and print the dependency graph "
        "of rules with their input and output files in the dot language. "
        "This is an intermediate solution between above DAG of jobs and the rule graph. "
        "Note that each rule is displayed once, hence the displayed graph will be "
        "cyclic if a rule appears in several steps of the workflow. "
        "Use this if above option leads to a DAG that is too large. "
        "Recommended use on Unix systems: snakemake --filegraph | dot | display. "
        "Note print statements in your Snakefile may interfere "
        "with visualization.",
    )
    group_utils.add_argument(
        "--d3dag",
        action="store_true",
        help="Print the DAG in D3.js compatible JSON format.",
    )
    group_utils.add_argument(
        "--summary",
        "-S",
        action="store_true",
        help="Print a summary of all files created by the workflow. The "
        "has the following columns: filename, modification time, "
        "rule version, status, plan.\n"
        "Thereby rule version contains the version "
        "the file was created with (see the version keyword of rules), and "
        "status denotes whether the file is missing, its input files are "
        "newer or if version or implementation of the rule changed since "
        "file creation. Finally the last column denotes whether the file "
        "will be updated or created during the next workflow execution.",
    )
    group_utils.add_argument(
        "--detailed-summary",
        "-D",
        action="store_true",
        help="Print a summary of all files created by the workflow. The "
        "has the following columns: filename, modification time, "
        "rule version, input file(s), shell command, status, plan.\n"
        "Thereby rule version contains the version "
        "the file was created with (see the version keyword of rules), and "
        "status denotes whether the file is missing, its input files are "
        "newer or if version or implementation of the rule changed since "
        "file creation. The input file and shell command columns are self "
        "explanatory. Finally the last column denotes whether the file "
        "will be updated or created during the next workflow execution.",
    )
    group_utils.add_argument(
        "--archive",
        metavar="FILE",
        type=Path,
        help="Archive the workflow into the given tar archive FILE. The archive "
        "will be created such that the workflow can be re-executed on a vanilla "
        "system. The function needs conda and git to be installed. "
        "It will archive every file that is under git version control. "
        "Note that it is best practice to have the Snakefile, config files, and "
        "scripts under version control. Hence, they will be included in the archive. "
        "Further, it will add input files that are not generated by "
        "by the workflow itself and conda environments. Note that symlinks are "
        "dereferenced. Supported formats are .tar, .tar.gz, .tar.bz2 and .tar.xz.",
    )
    group_utils.add_argument(
        "--cleanup-metadata",
        "--cm",
        nargs="+",
        metavar="FILE",
        type=Path,
        help="Cleanup the metadata "
        "of given files. That means that snakemake removes any tracked "
        "version info, and any marks that files are incomplete.",
    )
    group_utils.add_argument(
        "--cleanup-shadow",
        action="store_true",
        help="Cleanup old shadow directories which have not been deleted due "
        "to failures or power loss.",
    )
    group_utils.add_argument(
        "--skip-script-cleanup",
        action="store_true",
        help="Don't delete wrapper scripts used for execution",
    )
    group_utils.add_argument(
        "--unlock", action="store_true", help="Remove a lock on the working directory."
    )
    group_utils.add_argument(
        "--list-changes",
        "--lc",
        choices=ChangeType.all(),
        type=ChangeType.parse_choice,
        help="List all output files for which the given items (code, input, params) "
        "have changed since creation.",
    )
    group_utils.add_argument(
        "--list-input-changes",
        "--li",
        action="store_true",
        help="List all output files for which the defined input files have changed "
        "in the Snakefile (e.g. new input files were added in the rule "
        "definition or files were renamed). For listing input file "
        "modification in the filesystem, use `--summary`.",
    )
    group_utils.add_argument(
        "--list-params-changes",
        "--lp",
        action="store_true",
        help="List all output files for which the defined params have changed "
        "in the Snakefile.",
    )
    group_utils.add_argument(
        "--list-untracked",
        "--lu",
        action="store_true",
        help="List all files in the working directory that are not used in the  "
        "workflow. This can be used e.g. for identifying leftover files. Hidden files "
        "and directories are ignored.",
    )
    group_delete_output = group_utils.add_mutually_exclusive_group()
    group_delete_output.add_argument(
        "--delete-all-output",
        action="store_true",
        help="Remove all files generated by the workflow. Use together with `--dry-run` "
        "to list files without actually deleting anything. Note that this will "
        "not recurse into subworkflows. Write-protected files are not removed. "
        "Nevertheless, use with care!",
    )
    group_delete_output.add_argument(
        "--delete-temp-output",
        action="store_true",
        help="Remove all temporary files generated by the workflow. Use together "
        "with `--dry-run` to list files without actually deleting anything. Note "
        "that this will not recurse into subworkflows.",
    )
    group_utils.add_argument(
        "--keep-incomplete",
        action="store_true",
        help="Do not remove incomplete output files by failed jobs.",
    )
    group_utils.add_argument(
        "--drop-metadata",
        action="store_true",
        help="Drop metadata file tracking information after job finishes. "
        "Provenance-information based reports (e.g. `--report` and the "
        "`--list_x_changes` functions) will be empty or incomplete.",
    )
    group_utils.add_argument("--version", "-v", action="version", version=__version__)

    group_output = parser.add_argument_group("OUTPUT")
    group_output.add_argument(
        "--printshellcmds",
        "-p",
        action="store_true",
        help="Print out the shell commands that will be executed.",
    )
    group_output.add_argument(
        "--debug-dag",
        action="store_true",
        help="Print candidate and selected jobs (including their wildcards) while "
        "inferring DAG. This can help to debug unexpected DAG topology or errors.",
    )
    group_output.add_argument(
        "--nocolor", action="store_true", help="Do not use a colored output."
    )
    group_output.add_argument(
        "--quiet",
        "-q",
        nargs="*",
        choices=Quietness.choices(),
        default=None,
        parse_func=parse_quietness,
        help="Do not output certain information. "
        "If used without arguments, do not output any progress or rule "
        "information. Defining `all` results in no information being "
        "printed at all.",
    )
    group_output.add_argument(
        "--print-compilation",
        action="store_true",
        help="Print the python representation of the workflow.",
    )

    group_output.add_argument(
        "--verbose", action="store_true", help="Print debugging output."
    )

    group_behavior = parser.add_argument_group("BEHAVIOR")
    group_behavior.add_argument(
        "--force-use-threads",
        dest="force_use_threads",
        action="store_true",
        help="Force threads rather than processes. Helpful if shared memory (/dev/shm) is full or unavailable.",
    )
    group_behavior.add_argument(
        "--allow-ambiguity",
        "-a",
        action="store_true",
        help=(
            "Don't check for ambiguous rules and simply use the first if "
            "several can produce the same file. This allows the user to "
            "prioritize rules by their order in the snakefile."
        ),
    )
    group_behavior.add_argument(
        "--nolock", action="store_true", help="Do not lock the working directory"
    )
    group_behavior.add_argument(
        "--ignore-incomplete",
        "--ii",
        action="store_true",
        help="Do not check for incomplete output files.",
    )
    group_behavior.add_argument(
        "--max-inventory-time",
        type=int,
        default=20,
        metavar="SECONDS",
        help="Spend at most SECONDS seconds to create a file inventory for the working directory. "
        "The inventory vastly speeds up file modification and existence checks when computing "
        "which jobs need to be executed. However, creating the inventory itself can be slow, e.g. on "
        "network file systems. Hence, we do not spend more than a given amount of time and fall back "
        "to individual checks for the rest.",
    )
    group_behavior.add_argument(
        "--trust-io-cache",
        action="store_true",
        help=(
            "Tell Snakemake to assume that all input and output file existence and modification "
            "time queries performed in previous dryruns are still valid and therefore don't have to "
            "be repeated. This can lead to speed-ups, but implies that input and output have not "
            "been modified manually in between. Non dry-run execution will automatically "
            "invalidate the cache and lead to redoing the queries."
        ),
    )
    group_behavior.add_argument(
        "--max-checksum-file-size",
        default=1000000,
        metavar="SIZE",
        parse_func=parse_size_in_bytes,
        help=(
            "Compute the checksum during DAG computation and job postprocessing "
            "only for files that are smaller than the provided threshold (given in any valid size "
            "unit, e.g. 1MB, which is also the default). "
        ),
    )
    group_behavior.add_argument(
        "--latency-wait",
        "--output-wait",
        "-w",
        type=int,
        default=5,
        metavar="SECONDS",
        help="Wait given seconds if an output file of a job is not present after "
        "the job finished. This helps if your filesystem suffers from latency.",
    )
    group_behavior.add_argument(
        "--wait-for-free-local-storage",
        parse_func=parse_timespan,
        help=(
            "Wait for given timespan for enough free local storage when downloading "
            "remote storage files. If not set, no waiting is performed."
        ),
    )
    group_behavior.add_argument(
        "--wait-for-files",
        nargs="*",
        metavar="FILE",
        parse_func=set,
        help="Wait `--latency-wait` seconds for these "
        "files to be present before executing the workflow. "
        "This option is used internally to handle filesystem latency in cluster "
        "environments.",
    )
    group_behavior.add_argument(
        "--wait-for-files-file",
        metavar="FILE",
        help="Same behaviour as `--wait-for-files`, but file list is "
        "stored in file instead of being passed on the commandline. "
        "This is useful when the list of files is too long to be "
        "passed on the commandline.",
    )
    group_behavior.add_argument(
        "--queue-input-wait-time",
        metavar="SECONDS",
        type=int,
        default=10,
        help="Set the interval in seconds to check for new input in rules that use from_queue to obtain input files.",
    )
    group_behavior.add_argument(
        "--notemp",
        "--nt",
        action="store_true",
        help="Ignore temp() declarations. This is useful when running only "
        "a part of the workflow, since temp() would lead to deletion of "
        "probably needed files by other parts of the workflow.",
    )
    group_behavior.add_argument(
        "--all-temp",
        action="store_true",
        help="Mark all output files as temp files. This can be useful for CI testing, "
        "in order to save space.",
    )
    group_behavior.add_argument(
        "--unneeded-temp-files",
        parse_func=set,
        default=frozenset(),
        metavar="FILE",
        nargs="+",
        help="Given files will not be uploaded to storage and immediately deleted "
        "after job or group job completion.",
    )
    group_behavior.add_argument(
        "--keep-storage-local-copies",
        action="store_true",
        help="Keep local copies of remote input and output files.",
    )
    group_behavior.add_argument(
        "--not-retrieve-storage",
        action="store_true",
        help="Do not retrieve remote files (default is to retrieve remote files).",
    )
    group_behavior.add_argument(
        "--target-files-omit-workdir-adjustment",
        action="store_true",
        help="Do not adjust the paths of given target files relative to the working directory.",
    )
    group_behavior.add_argument(
        "--allowed-rules",
        nargs="+",
        help="Only consider given rules. If omitted, all rules in Snakefile are "
        "used. Note that this is intended primarily for internal use and may "
        "lead to unexpected results otherwise.",
    )
    group_behavior.add_argument(
        "--max-jobs-per-timespan",
        default="100/1s",
        type=MaxJobsPerTimespan.parse_choice,
        help="Maximal number of job submissions/executions per timespan. Format: <number><timespan>, e.g. 50/1m or 0.5/1s.",
    )
    group_behavior.add_argument(
        "--max-jobs-per-second",
        type=int,
        help="Maximal number of job submissions/executions per second. "
        "Deprecated in favor of `--max-jobs-per-timespan`.",
    )
    group_behavior.add_argument(
        "--max-status-checks-per-second",
        default=10,
        type=float,
        help="Maximal number of job status checks per second; fractions allowed.",
    )
    group_behavior.add_argument(
        "--seconds-between-status-checks",
        default=10,
        type=int,
        help="Number of seconds to wait between two rounds of status checks.",
    )
    group_behavior.add_argument(
        "--retries",
        "--restart-times",
        "-T",
        default=0,
        type=int,
        help="Number of times to restart failing jobs.",
    )
    group_behavior.add_argument(
        "--wrapper-prefix",
        default="https://github.com/snakemake/snakemake-wrappers/raw/",
        help="URL prefix for wrapper directive. Set this to use your fork or a local clone of the repository, "
        "e.g., use a git URL like `git+file://path/to/your/local/clone@`.",
    )
    group_behavior.add_argument(
        "--default-storage-provider",
        type=optional_str,
        help="Specify default storage provider to be used for "
        "all input and output files that don't yet specify "
        "one (e.g. `s3`). See https://snakemake.github.io/snakemake-plugin-catalog "
        "for available storage provider plugins. If not set or explicitly `none`, no "
        "default storage provider will be used.",
    )
    group_behavior.add_argument(
        "--default-storage-prefix",
        default="",
        help="Specify prefix for default storage provider. E.g. a bucket name.",
    )
    group_behavior.add_argument(
        "--local-storage-prefix",
        default=".snakemake/storage",
        type=maybe_base64(expandvars(Path)),
        help="Specify prefix for storing local copies of storage files and folders "
        "(e.g. local scratch disk). Environment variables will be expanded.",
    )
    group_behavior.add_argument(
        "--remote-job-local-storage-prefix",
        default=".snakemake/storage",
        type=maybe_base64(Path),
        help="Specify prefix for storing local copies of storage files and folders (e.g. local scratch disk) in "
        "case of remote jobs (e.g. cluster or cloud jobs). Environment variables will be expanded within the remote job.",
    )
    group_behavior.add_argument(
        "--shared-fs-usage",
        nargs="+",
        default=SharedFSUsage.all(),
        choices=SharedFSUsage.choices(),
        parse_func=SharedFSUsage.parse_choices_set,
        help="Set assumptions on shared filesystem for non-local "
        "workflow execution. To disable any sharing via the filesystem, "
        "specify `none`. "
        "Usually, the executor plugin sets this to a correct "
        "default. However, sometimes it is worth tuning this setting, e.g. for "
        "optimizing cluster performance. For example, when using "
        "`--default-storage-provider fs` together with a cluster executor like "
        "slurm, you might want to set "
        "`--shared-fs-usage persistence software-deployment sources source-cache`, "
        "such that software deployment and data provenance will be handled by NFS "
        "but input and output files will be handled exclusively by the storage "
        "provider.",
    )
    group_behavior.add_argument(
        "--scheduler-greediness",
        "--greediness",
        type=float,
        default=1.0,
        help="Set the greediness of scheduling. This value between 0 and 1 "
        "determines how careful jobs are selected for execution. The default "
        "value (1.0) provides the best speed and still acceptable scheduling "
        "quality.",
    )
    group_behavior.add_argument(
        "--scheduler-subsample",
        type=int,
        default=None,
        help="Set the number of jobs to be considered for scheduling. If number "
        "of ready jobs is greater than this value, this number of jobs is randomly "
        "chosen for scheduling; if number of ready jobs is lower, this option has "
        "no effect. This can be useful on very large DAGs, where the scheduler can "
        "take some time selecting which jobs to run.",
    )
    group_behavior.add_argument(
        "--no-hooks",
        action="store_true",
        help="Do not invoke onstart, onsuccess or onerror hooks after execution.",
    )
    group_behavior.add_argument(
        "--debug",
        action="store_true",
        help="Allow to debug rules with e.g. PDB. This flag "
        "allows to set breakpoints in run blocks.",
    )
    group_behavior.add_argument(
        "--runtime-profile",
        metavar="FILE",
        help="Profile Snakemake and write the output to FILE. This requires yappi "
        "to be installed.",
    )
    group_behavior.add_argument(
        "--local-groupid",
        default="local",
        help="Internal use only: Name for local groupid.",
    )
    group_behavior.add_argument(
        "--attempt",
        default=1,
        type=int,
        help="Internal use only: define the initial value of the attempt parameter.",
    )
    group_behavior.add_argument(
        "--show-failed-logs",
        action="store_true",
        help="Automatically display logs of failed jobs.",
    )

    group_behavior.add_argument(
        "--logger",
        nargs="+",
        default=[],
        choices=LoggerPluginRegistry().plugins.keys(),
        help="Specify one or more custom loggers, available via logger plugins.",
    )
    group_behavior.add_argument(
        "--job-deploy-sources",
        action="store_true",
        help="Whether the workflow sources shall be deployed before a remote job is "
        "started. Only applies if `--no-shared-fs` is set or executors are used that "
        "imply no shared FS (e.g. the kubernetes executor).",
    )
    group_behavior.add_argument(
        "--benchmark-extended",
        default=False,
        action="store_true",
        help="Write extended benchmarking metrics.",
    )

    group_cluster = parser.add_argument_group("REMOTE EXECUTION")

    group_cluster.add_argument(
        "--container-image",
        metavar="IMAGE",
        default=get_container_image(),
        help="Docker image to use, e.g., when submitting jobs to kubernetes. "
        "Defaults to https://hub.docker.com/r/snakemake/snakemake, tagged with "
        "the same version as the currently running Snakemake instance. "
        "Note that overwriting this value is up to your responsibility. "
        "Any used image has to contain a working snakemake installation "
        "that is compatible with (or ideally the same as) the currently "
        "running version.",
    )
    group_cluster.add_argument(
        "--immediate-submit",
        "--is",
        action="store_true",
        help="Immediately submit all jobs to the cluster instead of waiting "
        "for present input files. This will fail, unless you make "
        "the cluster aware of job dependencies, e.g. via: "
        "`$ snakemake --cluster 'sbatch --dependency {dependencies}'`. "
        "Assuming that your submit script (here sbatch) outputs the "
        "generated job id to the first stdout line, {dependencies} will "
        "be filled with space separated job ids this job depends on. "
        "Does not work for workflows that contain checkpoint rules.",
    )
    group_cluster.add_argument(
        "--jobscript",
        "--js",
        metavar="SCRIPT",
        help="Provide a custom job script for submission to the cluster. "
        "The default script resides as `jobscript.sh` in the installation directory.",
    )
    group_cluster.add_argument(
        "--jobname",
        "--jn",
        default="snakejob.{name}.{jobid}.sh",
        metavar="NAME",
        help="Provide a custom name for the jobscript that is submitted to the cluster (see `--cluster`). The wildcard `{jobid}` has to be present in the name.",
    )

    group_flux = parser.add_argument_group("FLUX")

    group_flux.add_argument(
        "--flux",
        action="store_true",
        help="Execute your workflow on a flux cluster. "
        "Flux can work with both a shared network filesystem (like NFS) or without. "
        "If you don't have a shared filesystem, additionally specify `--no-shared-fs`.",
    )

    group_deployment = parser.add_argument_group("SOFTWARE DEPLOYMENT")
    group_deployment.add_argument(
        "--software-deployment-method",
        "--deployment-method",
        "--deployment",
        "--sdm",
        nargs="+",
        choices=DeploymentMethod.choices(),
        parse_func=DeploymentMethod.parse_choices_set,
        default=set(),
        help="Specify software environment deployment method.",
    )
    group_deployment.add_argument(
        "--container-cleanup-images",
        action="store_true",
        help="Remove unused containers",
    )

    group_conda = parser.add_argument_group("CONDA")

    group_conda.add_argument(
        "--use-conda",
        action="store_true",
        help="If defined in the rule, run job in a conda environment. "
        "If this flag is not set, the conda directive is ignored.",
    )
    group_conda.add_argument(
        "--conda-not-block-search-path-envvars",
        action="store_true",
        help="Do not block environment variables that modify the search path "
        "(R_LIBS, PYTHONPATH, PERL5LIB, PERLLIB) when using conda environments.",
    )
    group_conda.add_argument(
        "--list-conda-envs",
        action="store_true",
        help="List all conda environments and their location on disk.",
    )
    group_conda.add_argument(
        "--conda-prefix",
        metavar="DIR",
        default=os.environ.get("SNAKEMAKE_CONDA_PREFIX", None),
        help="Specify a directory in which the `conda` and `conda-archive` "
        "directories are created. These are used to store conda environments "
        "and their archives, respectively. If not supplied, the value is set "
        "to the `.snakemake` directory relative to the invocation directory. "
        "If supplied, the `--use-conda` flag must also be set. The value may "
        "be given as a relative path, which will be extrapolated to the "
        "invocation directory, or as an absolute path. The value can also be "
        "provided via the environment variable $SNAKEMAKE_CONDA_PREFIX. "
        "In any case, the prefix may contain environment "
        "variables which will be properly expanded. "
        "Note that if you use remote execution "
        "e.g. on a cluster and you have node specific values for this, you should "
        "disable assuming shared fs for software-deployment (see `--shared-fs-usage`).",
    )
    group_conda.add_argument(
        "--conda-cleanup-envs",
        action="store_true",
        help="Cleanup unused conda environments.",
    )

    from snakemake.deployment.conda import CondaCleanupMode

    group_conda.add_argument(
        "--conda-cleanup-pkgs",
        type=CondaCleanupMode,
        const=CondaCleanupMode.tarballs,
        choices=list(CondaCleanupMode),
        default="tarballs",
        nargs="?",
        help="Cleanup conda packages after creating environments. "
        "In case of `tarballs` mode, will clean up all downloaded package tarballs. "
        "In case of `cache` mode, will additionally clean up unused package caches.",
    )
    group_conda.add_argument(
        "--conda-create-envs-only",
        action="store_true",
        help="If specified, only creates the job-specific "
        "conda environments then exits. The `--use-conda` "
        "flag must also be set.",
    )
    group_conda.add_argument(
        "--conda-frontend",
        default="conda",
        choices=["conda", "mamba"],
        help="Choose the conda frontend for installing environments.",
    )

    group_singularity = parser.add_argument_group("APPTAINER/SINGULARITY")

    group_singularity.add_argument(
        "--use-apptainer",
        "--use-singularity",
        action="store_true",
        help="If defined in the rule, run job within a apptainer/singularity container. "
        "If this flag is not set, the singularity directive is ignored.",
    )
    group_singularity.add_argument(
        "--apptainer-prefix",
        "--singularity-prefix",
        metavar="DIR",
        help="Specify a directory in which apptainer/singularity images will be stored."
        "If not supplied, the value is set "
        "to the `.snakemake` directory relative to the invocation directory. "
        "If supplied, the `--use-apptainer` flag must also be set. The value "
        "may be given as a relative path, which will be extrapolated to the "
        "invocation directory, or as an absolute path. If not supplied, "
        "APPTAINER_CACHEDIR is used. In any case, the prefix may contain environment "
        "variables which will be properly expanded. Note that if you use remote execution "
        "e.g. on a cluster and you have node specific values for this, you should "
        "disable assuming shared fs for software-deployment (see `--shared-fs-usage`).",
    )
    group_singularity.add_argument(
        "--apptainer-args",
        "--singularity-args",
        default="",
        metavar="ARGS",
        parse_func=maybe_base64(str),
        help="Pass additional args to apptainer/singularity.",
    )

    group_env_modules = parser.add_argument_group("ENVIRONMENT MODULES")

    group_env_modules.add_argument(
        "--use-envmodules",
        action="store_true",
        help="If defined in the rule, run job within the given environment "
        "modules, loaded in the given order. This can be combined with "
        "`--use-conda` and `--use-singularity`, which will then be only used as a "
        "fallback for rules which don't define environment modules.",
    )

    def help_internal(text):
        return f"Internal use only: {text}"

    group_internal = parser.add_argument_group("INTERNAL")
    group_internal.add_argument(
        "--scheduler-solver-path",
        help=help_internal("Set the PATH to search for scheduler solver binaries."),
    )
    group_internal.add_argument(
        "--deploy-sources",
        nargs=2,
        metavar=("QUERY", "CHECKSUM"),
        help=help_internal(
            "Deploy sources archive from given storage provider query to the current "
            "working subdirectory and control for archive checksum to proceed."
        ),
    )
    group_internal.add_argument(
        "--target-jobs",
        nargs="+",
        parse_func=parse_target_jobs_cli_args,
        default=set(),
        help=help_internal(
            "Target particular jobs by RULE:WILDCARD1=VALUE,WILDCARD2=VALUE,... "
        ),
    )
    group_internal.add_argument(
        "--mode",
        choices=ExecMode.all(),
        default=ExecMode.DEFAULT,
        type=ExecMode.parse_choice,
        help=help_internal("Set execution mode of Snakemake."),
    )

    # Add namespaced arguments to parser for each plugin
    ExecutorPluginRegistry().register_cli_args(parser)
    StoragePluginRegistry().register_cli_args(parser)
    ReportPluginRegistry().register_cli_args(parser)
    LoggerPluginRegistry().register_cli_args(parser)
    return parser


def generate_parser_metadata(parser, args):
    """Given a populated parser, generate the original command along with
    metadata that can be handed to a logger to use as needed.
    """
    command = "snakemake %s" % " ".join(
        parser._source_to_settings["command_line"][""][1]
    )
    metadata = args.__dict__
    metadata.update({"command": command})
    return metadata


def parse_args(argv):
    parser = get_argument_parser()
    args = parser.parse_args(argv)

    snakefile = resolve_snakefile(args.snakefile, allow_missing=True)
    workflow_profile = None
    if args.workflow_profile != "none":
        if args.workflow_profile:
            workflow_profile = args.workflow_profile
        elif snakefile is not None:
            # checking for default profile
            default_path = Path("profiles/default")
            workflow_profile_candidates = [
                default_path,
                Path(snakefile).parent.joinpath(default_path),
            ]
            for profile in workflow_profile_candidates:
                if profile.exists():
                    workflow_profile = profile
                    break

    if args.profile == "none":
        args.profile = None

    if (args.profile or workflow_profile) and args.mode == ExecMode.DEFAULT:
        # Reparse args while inferring config file from profile.
        # But only do this if the user has invoked Snakemake (ExecMode.DEFAULT)
        profiles = []
        if args.profile:
            profiles.append(args.profile)
        if workflow_profile:
            workflow_profile_stmt = f" {'and ' if profiles else ''}workflow specific profile {workflow_profile}"
            profiles.append(workflow_profile)
        else:
            workflow_profile_stmt = ""

        fmt_profiles = profiles if not workflow_profile else profiles[:-1]
        profile_stmt = (
            f"Using profile{'s' if len(profiles) > 1 else ''} {' and '.join(map(str, fmt_profiles))}"
            if fmt_profiles
            else "Using"
        )

        if profiles:
            print(
                f"{profile_stmt}{workflow_profile_stmt} for setting default command line arguments.",
                file=sys.stderr,
            )

        parser = get_argument_parser(profiles=profiles)
        args = parser.parse_args(argv)

    return parser, args


def parse_quietness(quietness) -> Set[Quietness]:
    if quietness is not None and len(quietness) == 0:
        # default case, set quiet to progress and rule
        quietness = {Quietness.PROGRESS, Quietness.RULES}
    else:
        quietness = Quietness.parse_choices_set(quietness)
    return quietness


def parse_edit_notebook(args):
    edit_notebook = None
    if args.draft_notebook:
        args.targets = {args.draft_notebook}
        edit_notebook = NotebookEditMode(draft_only=True)
    elif args.edit_notebook:
        args.targets = {args.edit_notebook}
        args.force = True
        edit_notebook = NotebookEditMode(args.notebook_listen)
    return edit_notebook


def parse_wait_for_files(args):
    from snakemake.io import wait_for_files

    aggregated_wait_for_files = args.wait_for_files
    if args.wait_for_files_file is not None:
        async_run(
            wait_for_files([args.wait_for_files_file], latency_wait=args.latency_wait)
        )

        with open(args.wait_for_files_file) as fd:
            extra_wait_files = [line.strip() for line in fd.readlines()]

        if aggregated_wait_for_files is None:
            aggregated_wait_for_files = extra_wait_files
        else:
            aggregated_wait_for_files.update(extra_wait_files)
    return aggregated_wait_for_files


def parse_rerun_triggers(values):
    return {RerunTrigger[x] for x in values}


def args_to_api(args, parser):
    """Convert argparse args to API calls."""

    # handle legacy executor names
    if args.dryrun:
        args.executor = "dryrun"
    elif args.touch:
        args.executor = "touch"
    elif args.executor is None:
        args.executor = "local"

    if args.report:
        args.reporter = "html"
        args.report_html_path = args.report
        args.report_html_stylesheet_path = args.report_stylesheet

    if args.report_after_run and args.report is None:
        raise CliException(
            "The option --report-after-run requires the --report option."
        )

    executor_plugin = ExecutorPluginRegistry().get_plugin(args.executor)
    executor_settings = executor_plugin.get_settings(args)

    storage_provider_settings = {
        name: StoragePluginRegistry().get_plugin(name).get_settings(args)
        for name in StoragePluginRegistry().get_registered_plugins()
    }

    log_handler_settings = {
        name: LoggerPluginRegistry().get_plugin(name).get_settings(args)
        for name in args.logger
    }

    if args.reporter:
        report_plugin = ReportPluginRegistry().get_plugin(args.reporter)
        report_settings = report_plugin.get_settings(args)
    else:
        report_plugin = None
        report_settings = None

    if args.cores is None:
        if executor_plugin.common_settings.local_exec:
            # use --jobs as an alias for --cores
            args.cores = args.jobs
            args.jobs = None
        elif executor_plugin.common_settings.dryrun_exec:
            args.cores = 1
            args.jobs = None

    # start profiler if requested
    if args.runtime_profile:
        import yappi

        yappi.start()

    edit_notebook = parse_edit_notebook(args)

    wait_for_files = parse_wait_for_files(args)

    with SnakemakeApi(
        OutputSettings(
            dryrun=args.dryrun,
            printshellcmds=args.printshellcmds,
            nocolor=args.nocolor,
            quiet=args.quiet,
            debug_dag=args.debug_dag,
            verbose=args.verbose,
            show_failed_logs=args.show_failed_logs,
            log_handler_settings=log_handler_settings,
            keep_logger=False,
            stdout=args.dryrun,
            benchmark_extended=args.benchmark_extended,
        )
    ) as snakemake_api:
        deployment_method = args.software_deployment_method
        if args.use_conda:
            deployment_method.add(DeploymentMethod.CONDA)
        if args.use_apptainer:
            deployment_method.add(DeploymentMethod.APPTAINER)
        if args.use_envmodules:
            deployment_method.add(DeploymentMethod.ENV_MODULES)

        try:
            storage_settings = StorageSettings(
                default_storage_provider=args.default_storage_provider,
                default_storage_prefix=args.default_storage_prefix,
                local_storage_prefix=args.local_storage_prefix,
                remote_job_local_storage_prefix=args.remote_job_local_storage_prefix,
                shared_fs_usage=args.shared_fs_usage,
                keep_storage_local=args.keep_storage_local_copies,
                retrieve_storage=not args.not_retrieve_storage,
                notemp=args.notemp,
                all_temp=args.all_temp,
                unneeded_temp_files=args.unneeded_temp_files,
                wait_for_free_local_storage=args.wait_for_free_local_storage,
            )

            if args.deploy_sources:
                query, checksum = args.deploy_sources
                snakemake_api.deploy_sources(
                    query,
                    checksum,
                    storage_settings=storage_settings,
                    storage_provider_settings=storage_provider_settings,
                )
            else:
                workflow_api = snakemake_api.workflow(
                    resource_settings=ResourceSettings(
                        cores=args.cores,
                        nodes=args.jobs,
                        local_cores=args.local_cores,
                        max_threads=args.max_threads,
                        resources=args.resources,
                        overwrite_threads=args.set_threads,
                        overwrite_scatter=args.set_scatter,
                        overwrite_resource_scopes=args.set_resource_scopes,
                        overwrite_resources=args.set_resources,
                        default_resources=args.default_resources,
                    ),
                    config_settings=ConfigSettings(
                        config=parse_config(args.config),
                        configfiles=args.configfile,
                        config_args=args.config,
                        replace_workflow_config=args.replace_workflow_config,
                    ),
                    storage_settings=storage_settings,
                    storage_provider_settings=storage_provider_settings,
                    workflow_settings=WorkflowSettings(
                        wrapper_prefix=args.wrapper_prefix,
                        exec_mode=args.mode,
                        cache=args.cache,
                        consider_ancient=args.consider_ancient,
                    ),
                    deployment_settings=DeploymentSettings(
                        deployment_method=deployment_method,
                        conda_prefix=args.conda_prefix,
                        conda_cleanup_pkgs=args.conda_cleanup_pkgs,
                        conda_base_path=args.conda_base_path,
                        conda_frontend=args.conda_frontend,
                        conda_not_block_search_path_envvars=args.conda_not_block_search_path_envvars,
                        apptainer_args=args.apptainer_args,
                        apptainer_prefix=args.apptainer_prefix,
                    ),
                    snakefile=args.snakefile,
                    workdir=args.directory,
                )

                if args.lint:
                    any_lint = workflow_api.lint()
                    if any_lint:
                        # trigger exit code 1
                        return False
                elif args.list_target_rules:
                    workflow_api.list_rules(only_targets=True)
                elif args.list_rules:
                    workflow_api.list_rules(only_targets=False)
                elif args.print_compilation:
                    workflow_api.print_compilation()
                else:

                    print_dag_as = None
                    if args.dag:
                        print_dag_as = args.dag
                    elif args.rulegraph:
                        print_dag_as = args.rulegraph

                    dag_api = workflow_api.dag(
                        dag_settings=DAGSettings(
                            targets=args.targets,
                            target_jobs=args.target_jobs,
                            target_files_omit_workdir_adjustment=args.target_files_omit_workdir_adjustment,
                            batch=args.batch,
                            forcetargets=args.force,
                            forceall=args.forceall,
                            forcerun=args.forcerun,
                            until=args.until,
                            omit_from=args.omit_from,
                            force_incomplete=args.rerun_incomplete,
                            allowed_rules=args.allowed_rules,
                            rerun_triggers=args.rerun_triggers,
                            max_inventory_wait_time=args.max_inventory_time,
                            trust_io_cache=args.trust_io_cache,
                            max_checksum_file_size=args.max_checksum_file_size,
                            strict_evaluation=args.strict_dag_evaluation,
                            print_dag_as=print_dag_as,
                        ),
                    )

                    # None indicates that the flag is not being used, a set indicates used
                    if args.preemptible_rules is not None:
                        # If no explicit rules provided (empty set) we assume all preemptible
                        # Checking for "not" evaluates for an empty set
                        if not args.preemptible_rules:
                            preemptible_rules = PreemptibleRules(all=True)
                        else:
                            preemptible_rules = PreemptibleRules(
                                rules=args.preemptible_rules
                            )
                    else:
                        preemptible_rules = PreemptibleRules()

                    if args.containerize:
                        dag_api.containerize()
                    elif report_plugin is not None and not args.report_after_run:
                        dag_api.create_report(
                            reporter=args.reporter,
                            report_settings=report_settings,
                        )
                    elif args.generate_unit_tests:
                        dag_api.generate_unit_tests(args.generate_unit_tests)
                    elif args.dag:
                        dag_api.printdag()
                    elif args.rulegraph:
                        dag_api.printrulegraph()
                    elif args.filegraph:
                        dag_api.printfilegraph()
                    elif args.d3dag:
                        dag_api.printd3dag()
                    elif args.unlock:
                        dag_api.unlock()
                    elif args.cleanup_metadata:
                        dag_api.cleanup_metadata(args.cleanup_metadata)
                    elif args.conda_cleanup_envs:
                        dag_api.conda_cleanup_envs()
                    elif args.conda_create_envs_only:
                        dag_api.conda_create_envs()
                    elif args.list_conda_envs:
                        dag_api.conda_list_envs()
                    elif args.cleanup_shadow:
                        dag_api.cleanup_shadow()
                    elif args.container_cleanup_images:
                        dag_api.container_cleanup_images()
                    elif args.list_changes:
                        dag_api.list_changes(args.list_changes)
                    elif args.list_input_changes:
                        dag_api.list_changes(ChangeType.INPUT)
                    elif args.list_params_changes:
                        dag_api.list_changes(ChangeType.PARAMS)
                    elif args.list_untracked:
                        dag_api.list_untracked()
                    elif args.summary:
                        dag_api.summary()
                    elif args.detailed_summary:
                        dag_api.summary(detailed=True)
                    elif args.archive:
                        dag_api.archive(args.archive)
                    elif args.delete_all_output:
                        dag_api.delete_output(dryrun=args.dryrun)
                    elif args.delete_temp_output:
                        dag_api.delete_output(only_temp=True, dryrun=args.dryrun)
                    else:
                        dag_api.execute_workflow(
                            executor=args.executor,
                            execution_settings=ExecutionSettings(
                                keep_going=args.keep_going,
                                debug=args.debug,
                                standalone=True,
                                ignore_ambiguity=args.allow_ambiguity,
                                lock=not args.nolock,
                                ignore_incomplete=args.ignore_incomplete,
                                latency_wait=args.latency_wait,
                                wait_for_files=wait_for_files,
                                no_hooks=args.no_hooks,
                                retries=args.retries,
                                attempt=args.attempt,
                                use_threads=args.force_use_threads,
                                shadow_prefix=args.shadow_prefix,
                                keep_incomplete=args.keep_incomplete,
                                keep_metadata=not args.drop_metadata,
                                edit_notebook=edit_notebook,
                                cleanup_scripts=not args.skip_script_cleanup,
                                queue_input_wait_time=args.queue_input_wait_time,
                            ),
                            remote_execution_settings=RemoteExecutionSettings(
                                jobname=args.jobname,
                                jobscript=args.jobscript,
                                max_status_checks_per_second=args.max_status_checks_per_second,
                                seconds_between_status_checks=args.seconds_between_status_checks,
                                container_image=args.container_image,
                                preemptible_retries=args.preemptible_retries,
                                preemptible_rules=preemptible_rules,
                                envvars=args.envvars,
                                immediate_submit=args.immediate_submit,
                                job_deploy_sources=args.job_deploy_sources,
                                precommand=args.precommand,
                            ),
                            scheduling_settings=SchedulingSettings(
                                prioritytargets=args.prioritize,
                                scheduler=args.scheduler,
                                ilp_solver=args.scheduler_ilp_solver,
                                solver_path=args.scheduler_solver_path,
                                greediness=args.scheduler_greediness,
                                subsample=args.scheduler_subsample,
                                max_jobs_per_second=args.max_jobs_per_second,
                                max_jobs_per_timespan=args.max_jobs_per_timespan,
                            ),
                            group_settings=GroupSettings(
                                group_components=args.group_components,
                                overwrite_groups=args.groups,
                                local_groupid=args.local_groupid,
                            ),
                            executor_settings=executor_settings,
                        )

                        if report_plugin is not None and args.report_after_run:
                            dag_api.create_report(
                                reporter=args.reporter,
                                report_settings=report_settings,
                            )

        except Exception as e:
            snakemake_api.print_exception(e)
            return False

    # store profiler results if requested
    if args.runtime_profile:
        with open(args.runtime_profile, "w") as out:
            profile = yappi.get_func_stats()
            profile.sort("totaltime")
            profile.print_all(
                out=out,
                columns={
                    0: ("name", 120),
                    1: ("ncall", 10),
                    2: ("tsub", 8),
                    3: ("ttot", 8),
                    4: ("tavg", 8),
                },
            )
    return True


def main(argv=None):
    """Main entry point."""

    try:
        parser, args = parse_args(argv)
        success = args_to_api(args, parser)
    except Exception as e:
        print_exception(e)
        sys.exit(1)
    sys.exit(0 if success else 1)
