__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import sys

from snakemake import logging
from snakemake.api import snakemake

if sys.version_info < (3, 9):
    raise ValueError("Snakemake requires at least Python 3.9.")

import os
import glob
from argparse import ArgumentDefaultsHelpFormatter
import logging as _logging
from pathlib import Path
import re
import threading
import webbrowser
from functools import partial
import importlib
import shlex
from importlib.machinery import SourceFileLoader

from snakemake_interface_executor_plugins.utils import url_can_parse, ExecMode
from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry

from snakemake.target_jobs import parse_target_jobs_cli_args
from snakemake.workflow import Workflow
from snakemake.dag import Batch
from snakemake.exceptions import (
    CliException,
    ResourceScopesException,
    print_exception,
    WorkflowError,
)
from snakemake.logging import setup_logger, logger, SlackLogger, WMSLogger
from snakemake.io import load_configfile, wait_for_files
from snakemake.shell import shell
from snakemake.utils import update_config, available_cpu_count
from snakemake.common import (
    RERUN_TRIGGERS,
    __version__,
    MIN_PY_VERSION,
    get_appdirs,
    dict_to_key_value_args,
    parse_key_value_arg,
)
from snakemake.resources import ResourceScopes, parse_resources, DefaultResources

SNAKEFILE_CHOICES = [
    "Snakefile",
    "snakefile",
    "workflow/Snakefile",
    "workflow/snakefile",
]



def parse_set_threads(args):
    return parse_set_ints(
        args.set_threads,
        "Invalid threads definition: entries have to be defined as RULE=THREADS pairs "
        "(with THREADS being a positive integer).",
    )


def parse_set_resources(args):
    errmsg = (
        "Invalid resource definition: entries have to be defined as RULE:RESOURCE=VALUE, with "
        "VALUE being a positive integer or a string."
    )

    from collections import defaultdict

    assignments = defaultdict(dict)
    if args.set_resources is not None:
        for entry in args.set_resources:
            key, value = parse_key_value_arg(entry, errmsg=errmsg)
            key = key.split(":")
            if not len(key) == 2:
                raise ValueError(errmsg)
            rule, resource = key
            try:
                value = int(value)
            except ValueError:
                assignments[rule][resource] = value
                continue
            if value < 0:
                raise ValueError(errmsg)
            assignments[rule][resource] = value
    return assignments


def parse_set_scatter(args):
    return parse_set_ints(
        args.set_scatter,
        "Invalid scatter definition: entries have to be defined as NAME=SCATTERITEMS pairs "
        "(with SCATTERITEMS being a positive integer).",
    )


def parse_set_resource_scope(args):
    err_msg = (
        "Invalid resource scopes: entries must be defined as RESOURCE=SCOPE pairs, "
        "where SCOPE is either 'local', 'global', or 'excluded'"
    )
    if args.set_resource_scopes is not None:
        try:
            return ResourceScopes(
                parse_key_value_arg(entry, errmsg=err_msg)
                for entry in args.set_resource_scopes
            )
        except ResourceScopesException as err:
            invalid_resources = ", ".join(
                f"'{res}={scope}'" for res, scope in err.invalid_resources.items()
            )
            raise ValueError(f"{err.msg} (got {invalid_resources})")

    return ResourceScopes()


def parse_set_ints(arg, errmsg):
    assignments = dict()
    if arg is not None:
        for entry in arg:
            key, value = parse_key_value_arg(entry, errmsg=errmsg)
            try:
                value = int(value)
            except ValueError:
                raise ValueError(errmsg)
            if value < 0:
                raise ValueError(errmsg)
            assignments[key] = value
    return assignments


def parse_batch(args):
    errmsg = "Invalid batch definition: batch entry has to be defined as RULE=BATCH/BATCHES (with integers BATCH <= BATCHES, BATCH >= 1)."
    if args.batch is not None:
        rule, batchdef = parse_key_value_arg(args.batch, errmsg=errmsg)
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
    if args.groups is not None:
        for entry in args.groups:
            rule, group = parse_key_value_arg(entry, errmsg=errmsg)
            overwrite_groups[rule] = group
    return overwrite_groups


def parse_group_components(args):
    errmsg = "Invalid group components definition: entries have to be defined as GROUP=COMPONENTS pairs (with COMPONENTS being a positive integer)"
    group_components = dict()
    if args.group_components is not None:
        for entry in args.group_components:
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


def parse_config(args):
    """Parse config from args."""
    import yaml

    yaml_base_load = lambda s: yaml.load(s, Loader=yaml.loader.BaseLoader)
    parsers = [int, float, _bool_parser, yaml_base_load, str]
    config = dict()
    if args.config is not None:
        valid = re.compile(r"[a-zA-Z_]\w*$")
        for entry in args.config:
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


def parse_cores(cores, allow_none=False):
    if cores is None:
        if allow_none:
            return cores
        raise CliException(
            "Error: you need to specify the maximum number of CPU cores to "
            "be used at the same time. If you want to use N cores, say --cores N "
            "or -cN. For all cores on your system (be sure that this is "
            "appropriate) use --cores all. For no parallelization use --cores 1 or "
            "-c1."
        )
    if cores == "all":
        return available_cpu_count()
    try:
        return int(cores)
    except ValueError:
        raise CliException(
            "Error parsing number of cores (--cores, -c, -j): must be integer, "
            "empty, or 'all'."
        )


def parse_jobs(jobs, allow_none=False):
    if jobs is None:
        if allow_none:
            return jobs
        raise CliException(
            "Error: you need to specify the maximum number of jobs to "
            "be queued or executed at the same time with --jobs or -j."
        )
    if jobs == "unlimited":
        return sys.maxsize
    try:
        return int(jobs)
    except ValueError:
        raise CliException(
            "Error parsing number of jobs (--jobs, -j): must be integer."
        )


def parse_cores_jobs(cores, jobs, no_exec, non_local_exec, dryrun):
    if no_exec or dryrun:
        cores = parse_cores(cores, allow_none=True) or 1
        jobs = parse_jobs(jobs, allow_none=True) or 1
    elif non_local_exec:
        cores = parse_cores(cores, allow_none=True)
        jobs = parse_jobs(jobs)
    else:
        cores = parse_cores(cores or jobs)
        jobs = None

    return cores, jobs


def get_profile_file(profile, file, return_default=False):
    dirs = get_appdirs()
    if os.path.exists(profile):
        search_dirs = [os.path.dirname(profile)]
        profile = os.path.basename(profile)
    else:
        search_dirs = [os.getcwd(), dirs.user_config_dir, dirs.site_config_dir]
    get_path = lambda d: os.path.join(d, profile, file)
    for d in search_dirs:
        p = get_path(d)
        # "file" can actually be a full command. If so, `p` won't exist as the
        # below would check if e.g. '/path/to/profile/script --arg1 val --arg2'
        # exists. To fix this, we use shlex.split() to get the path to the
        # script. We check for both, in case the path contains spaces or some
        # other thing that would cause shlex.split() to mangle the path
        # inaccurately.
        if os.path.exists(p) or os.path.exists(shlex.split(p)[0]):
            return p

    if return_default:
        return file
    return None


def get_argument_parser(profiles=None):
    """Generate and return argument parser."""
    import configargparse
    from snakemake.profiles import ProfileConfigFileParser

    dirs = get_appdirs()
    config_files = []
    if profiles:
        for profile in profiles:
            if profile == "":
                print("Error: invalid profile name.", file=sys.stderr)
                exit(1)

            config_file = get_profile_file(profile, "config.yaml")
            if config_file is None:
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
            config_files.append(config_file)

    parser = configargparse.ArgumentParser(
        description="Snakemake is a Python based language and execution "
        "environment for GNU Make-like workflows.",
        formatter_class=ArgumentDefaultsHelpFormatter,
        default_config_files=config_files,
        config_file_parser_class=ProfileConfigFileParser,
    )

    group_exec = parser.add_argument_group("EXECUTION")

    group_exec.add_argument(
        "target",
        nargs="*",
        default=None,
        help="Targets to build. May be rules or files.",
    )

    group_exec.add_argument(
        "--dry-run",
        "--dryrun",
        "-n",
        dest="dryrun",
        action="store_true",
        help="Do not execute anything, and display what would be done. "
        "If you have a very large workflow, use --dry-run --quiet to just "
        "print a summary of the DAG of jobs.",
    )

    group_exec.add_argument(
        "--profile",
        help=f"""
            Name of profile to use for configuring
            Snakemake. Snakemake will search for a corresponding
            folder in {dirs.site_config_dir} and {dirs.user_config_dir}. Alternatively, this can be an
            absolute or relative path.
            The profile folder has to contain a file 'config.yaml'.
            This file can be used to set default values for command
            line options in YAML format. For example,
            '--cluster qsub' becomes 'cluster: qsub' in the YAML
            file. Profiles can be obtained from
            https://github.com/snakemake-profiles.
            The profile can also be set via the environment variable $SNAKEMAKE_PROFILE.
            To override this variable and use no profile at all, provide the value 'none'
            to this argument.
            """,
        env_var="SNAKEMAKE_PROFILE",
    )

    group_exec.add_argument(
        "--workflow-profile",
        help="""
            Path (relative to current directory) to workflow specific profile 
            folder to use for configuring Snakemake with parameters specific for this
            workflow (like resources).
            If this flag is not used, Snakemake will by default use 
            'profiles/default' if present (searched both relative to current directory
            and relative to Snakefile, in this order).
            For skipping any workflow specific profile provide the special value 'none'.
            Settings made in the workflow profile will override settings made in the
            general profile (see --profile).
            The profile folder has to contain a file 'config.yaml'.
            This file can be used to set default values for command
            line options in YAML format. For example,
            '--cluster qsub' becomes 'cluster: qsub' in the YAML
            file. It is advisable to use the workflow profile to set
            or overwrite e.g. workflow specific resources like the amount of threads
            of a particular rule or the amount of memory needed.
            Note that in such cases, the arguments may be given as nested YAML mappings 
            in the profile, e.g. 'set-threads: myrule: 4' instead of 'set-threads: myrule=4'.
            """,
        env_var="SNAKEMAKE_PROFILE",
    )

    group_exec.add_argument(
        "--cache",
        nargs="*",
        metavar="RULE",
        help="Store output files of given rules in a central cache given by the environment "
        "variable $SNAKEMAKE_OUTPUT_CACHE. Likewise, retrieve output files of the given rules "
        "from this cache if they have been created before (by anybody writing to the same cache), "
        "instead of actually executing the rules. Output files are identified by hashing all "
        "steps, parameters and software stack (conda envs or containers) needed to create them.",
    )

    group_exec.add_argument(
        "--snakefile",
        "-s",
        metavar="FILE",
        help=(
            "The workflow definition in form of a snakefile."
            "Usually, you should not need to specify this. "
            "By default, Snakemake will search for {} "
            "beneath the current working "
            "directory, in this order. "
            "Only if you definitely want a different layout, "
            "you need to use this parameter."
        ).format(", ".join(map("'{}'".format, SNAKEFILE_CHOICES))),
    )
    group_exec.add_argument(
        "--cores",
        "-c",
        action="store",
        const=available_cpu_count(),
        nargs="?",
        metavar="N",
        help=(
            "Use at most N CPU cores/jobs in parallel. "
            "If N is omitted or 'all', the limit is set to the number of "
            "available CPU cores. "
            "In case of cluster/cloud execution, this argument sets the maximum number "
            "of cores requested from the cluster or cloud scheduler. (See "
            "https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#"
            "resources-remote-execution for more info)"
            "This number is available to rules via workflow.cores."
        ),
    )
    group_exec.add_argument(
        "--jobs",
        "-j",
        metavar="N",
        nargs="?",
        const=available_cpu_count(),
        action="store",
        help=(
            "Use at most N CPU cluster/cloud jobs in parallel. For local execution this is "
            "an alias for --cores. Note: Set to 'unlimited' in case, this does not play a role."
        ),
    )
    group_exec.add_argument(
        "--local-cores",
        action="store",
        default=available_cpu_count(),
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
        nargs="*",
        metavar="NAME=INT",
        help=(
            "Define additional resources that shall constrain the scheduling "
            "analogously to --cores (see above). A resource is defined as "
            "a name and an integer value. E.g. --resources mem_mb=1000. Rules can "
            "use resources by defining the resource keyword, e.g. "
            "resources: mem_mb=600. If now two rules require 600 of the resource "
            "'mem_mb' they won't be run in parallel by the scheduler. In "
            "cluster/cloud mode, this argument will also constrain the amount of "
            "resources requested from the server. (See "
            "https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#"
            "resources-remote-execution for more info)"
        ),
    )
    group_exec.add_argument(
        "--set-threads",
        metavar="RULE=THREADS",
        nargs="+",
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
        "overwriting rules individually with --set-threads.",
    )
    group_exec.add_argument(
        "--set-resources",
        metavar="RULE:RESOURCE=VALUE",
        nargs="+",
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
        help="Overwrite number of scatter items of scattergather processes. This allows to fine-tune "
        "workflow parallelization. Thereby, SCATTERITEMS has to be a positive integer, and NAME has to be "
        "the name of the scattergather process defined via a scattergather directive in the workflow.",
    )
    group_exec.add_argument(
        "--set-resource-scopes",
        metavar="RESOURCE=[global|local]",
        nargs="+",
        help="Overwrite resource scopes. A scope determines how a constraint is "
        "reckoned in cluster execution. With RESOURCE=local, a constraint applied to "
        "RESOURCE using --resources will be considered the limit for each group "
        "submission. With RESOURCE=global, the constraint will apply across all groups "
        "cumulatively. By default, only `mem_mb` and `disk_mb` are considered local, "
        "all other resources are global. This may be modified in the snakefile using "
        "the `resource_scopes:` directive. Note that number of threads, specified via "
        "--cores, is always considered local. (See "
        "https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#"
        "resources-remote-execution for more info)",
    )
    group_exec.add_argument(
        "--default-resources",
        "--default-res",
        nargs="*",
        metavar="NAME=INT",
        help=(
            "Define default values of resources for rules that do not define their own values. "
            "In addition to plain integers, python expressions over inputsize are allowed (e.g. '2*input.size_mb'). "
            "The inputsize is the sum of the sizes of all input files of a rule. "
            "By default, Snakemake assumes a default for mem_mb, disk_mb, and tmpdir (see below). "
            "This option allows to add further defaults (e.g. account and partition for slurm) or to overwrite these default values. "
            "The defaults are 'mem_mb=max(2*input.size_mb, 1000)', "
            "'disk_mb=max(2*input.size_mb, 1000)' "
            "(i.e., default disk and mem usage is twice the input file size but at least 1GB), and "
            "the system temporary directory (as given by $TMPDIR, $TEMP, or $TMP) is used for the tmpdir resource. "
            "The tmpdir resource is automatically used by shell commands, scripts and wrappers to store temporary data (as it is "
            "mirrored into $TMPDIR, $TEMP, and $TMP for the executed subprocesses). "
            "If this argument is not specified at all, Snakemake just uses the tmpdir resource as outlined above."
        ),
    )

    group_exec.add_argument(
        "--preemption-default",
        type=int,
        default=None,
        help=(
            "A preemptible instance can be requested when using the Google Life Sciences API. If you set a --preemption-default,"
            "all rules will be subject to the default. Specifically, this integer is the number of restart attempts that will be "
            "made given that the instance is killed unexpectedly. Note that preemptible instances have a maximum running time of 24 "
            "hours. If you want to set preemptible instances for only a subset of rules, use --preemptible-rules instead."
        ),
    )

    group_exec.add_argument(
        "--preemptible-rules",
        nargs="+",
        default=None,
        help=(
            "A preemptible instance can be requested when using the Google Life Sciences API. If you want to use these instances "
            "for a subset of your rules, you can use --preemptible-rules and then specify a list of rule and integer pairs, where "
            "each integer indicates the number of restarts to use for the rule's instance in the case that the instance is "
            "terminated unexpectedly. --preemptible-rules can be used in combination with --preemption-default, and will take "
            "priority. Note that preemptible instances have a maximum running time of 24. If you want to apply a consistent "
            "number of retries across all your rules, use --preemption-default instead. "
            "Example: snakemake --preemption-default 10 --preemptible-rules map_reads=3 call_variants=0"
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
            "the workflow. Default values can be set by providing a JSON file "
            "(see Documentation)."
        ),
    )
    group_exec.add_argument(
        "--configfile",
        "--configfiles",
        nargs="+",
        metavar="FILE",
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
        "--envvars",
        nargs="+",
        metavar="VARNAME",
        help="Environment variables to pass to cloud jobs.",
    )
    group_exec.add_argument(
        "--directory",
        "-d",
        metavar="DIR",
        action="store",
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
            "--force, --forceall, or --forcerun. Note however that you lose "
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
        choices=RERUN_TRIGGERS,
        default=RERUN_TRIGGERS,
        help="Define what triggers the rerunning of a job. By default, "
        "all triggers are used, which guarantees that results are "
        "consistent with the workflow code and configuration. If you "
        "rather prefer the traditional way of just considering "
        "file modification dates, use '--rerun-trigger mtime'.",
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
        choices=ExecutorPluginRegistry().plugins,
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
        help=(
            "Force the re-execution or creation of the given rules or files."
            " Use this option if you changed a rule and want to have all its "
            "output in your workflow updated."
        ),
    )
    group_exec.add_argument(
        "--prioritize",
        "-P",
        nargs="+",
        metavar="TARGET",
        help=(
            "Tell the scheduler to assign creation of given targets "
            "(and all their dependencies) highest priority. (EXPERIMENTAL)"
        ),
    )
    group_exec.add_argument(
        "--batch",
        metavar="RULE=BATCH/BATCHES",
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
            "Specify a directory in which the 'shadow' directory is created. "
            "If not supplied, the value is set to the '.snakemake' directory relative "
            "to the working directory."
        ),
    )

    try:
        import pulp

        lp_solvers = pulp.list_solvers(onlyAvailable=True)
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
        "--wms-monitor",
        action="store",
        nargs="?",
        help=(
            "IP and port of workflow management system to monitor the execution of snakemake (e.g. http://127.0.0.1:5000)"
            " Note that if your service requires an authorization token, you must export WMS_MONITOR_TOKEN in the environment."
        ),
    )
    group_exec.add_argument(
        "--wms-monitor-arg",
        nargs="*",
        metavar="NAME=VALUE",
        help=(
            "If the workflow management service accepts extra arguments, provide."
            " them in key value pairs with --wms-monitor-arg. For example, to run"
            " an existing workflow using a wms monitor, you can provide the pair "
            " id=12345 and the arguments will be provided to the endpoint to "
            " first interact with the workflow"
        ),
    )
    group_exec.add_argument(
        "--scheduler-ilp-solver",
        default=recommended_lp_solver,
        choices=lp_solvers,
        help=("Specifies solver to be utilized when selecting ilp-scheduler."),
    )
    group_exec.add_argument(
        "--scheduler-solver-path",
        help="Set the PATH to search for scheduler solver binaries (internal use only).",
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

    # TODO add group_partitioning, allowing to define --group rulename=groupname.
    # i.e. setting groups via the CLI for improving cluster performance given
    # available resources.
    # TODO add an additional flag --group-components groupname=3, allowing to set the
    # number of connected components a group is allowed to span. By default, this is 1
    # (as now), but the flag allows to extend this. This can be used to run e.g.
    # 3 jobs of the same rule in the same group, although they are not connected.
    # Can be helpful for putting together many small jobs or benefitting of shared memory
    # setups.

    group_group = parser.add_argument_group("GROUPING")
    group_group.add_argument(
        "--groups",
        nargs="+",
        help="Assign rules to groups (this overwrites any "
        "group definitions from the workflow).",
    )
    group_group.add_argument(
        "--group-components",
        nargs="+",
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
        help="Create an HTML report with results and statistics. "
        "This can be either a .html file or a .zip file. "
        "In the former case, all results are embedded into the .html (this only works for small data). "
        "In the latter case, results are stored along with a file report.html in the zip archive. "
        "If no filename is given, an embedded report.html is the default.",
    )
    group_report.add_argument(
        "--report-stylesheet",
        metavar="CSSFILE",
        help="Custom stylesheet to use for report. In particular, this can be used for "
        "branding the report with e.g. a custom logo, see docs.",
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
        "closing the notebook and hitting the 'Quit' button on the jupyter dashboard. "
        "Afterwards, the updated notebook will be automatically stored in the path defined in the rule. "
        "If the notebook is not yet present, this will create an empty draft. ",
    )
    group_notebooks.add_argument(
        "--notebook-listen",
        metavar="IP:PORT",
        default="localhost:8888",
        help="The IP address and PORT the notebook server used for editing the notebook (--edit-notebook) will listen on.",
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
        help="Automatically generate unit tests for each workflow rule. "
        "This assumes that all input files of each job are already present. "
        "Rules without a job with present input files will be skipped (a warning will be issued). "
        "For each rule, one test case will be "
        "created in the specified test folder (.tests/unit by default). After "
        "successful execution, tests can be run with "
        "'pytest TESTPATH'.",
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
    group_utils.add_argument(
        "--dag",
        action="store_true",
        help="Do not execute anything and print the directed "
        "acyclic graph of jobs in the dot language. Recommended "
        "use on Unix systems: snakemake --dag | dot | display. "
        "Note print statements in your Snakefile may interfere "
        "with visualization.",
    )
    group_utils.add_argument(
        "--rulegraph",
        action="store_true",
        help="Do not execute anything and print the dependency graph "
        "of rules in the dot language. This will be less "
        "crowded than above DAG of jobs, but also show less information. "
        "Note that each rule is displayed once, hence the displayed graph will be "
        "cyclic if a rule appears in several steps of the workflow. "
        "Use this if above option leads to a DAG that is too large. "
        "Recommended use on Unix systems: snakemake --rulegraph | dot | display. "
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
        "Thereby rule version contains the version"
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
        help="Archive the workflow into the given tar archive FILE. The archive "
        "will be created such that the workflow can be re-executed on a vanilla "
        "system. The function needs conda and git to be installed. "
        "It will archive every file that is under git version control. "
        "Note that it is best practice to have the Snakefile, config files, and "
        "scripts under version control. Hence, they will be included in the archive. "
        "Further, it will add input files that are not generated by "
        "by the workflow itself and conda environments. Note that symlinks are "
        "dereferenced. Supported "
        "formats are .tar, .tar.gz, .tar.bz2 and .tar.xz.",
    )
    group_utils.add_argument(
        "--cleanup-metadata",
        "--cm",
        nargs="+",
        metavar="FILE",
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
        "--list-version-changes",
        "--lv",
        action="store_true",
        help="List all output files that have been created with "
        "a different version (as determined by the version keyword).",
    )
    group_utils.add_argument(
        "--list-code-changes",
        "--lc",
        action="store_true",
        help="List all output files for which the rule body (run or shell) have "
        "changed in the Snakefile.",
    )
    group_utils.add_argument(
        "--list-input-changes",
        "--li",
        action="store_true",
        help="List all output files for which the defined input files have changed "
        "in the Snakefile (e.g. new input files were added in the rule "
        "definition or files were renamed). For listing input file "
        "modification in the filesystem, use --summary.",
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
    group_utils.add_argument(
        "--delete-all-output",
        action="store_true",
        help="Remove all files generated by the workflow. Use together with --dry-run "
        "to list files without actually deleting anything. Note that this will "
        "not recurse into subworkflows. Write-protected files are not removed. "
        "Nevertheless, use with care!",
    )
    group_utils.add_argument(
        "--delete-temp-output",
        action="store_true",
        help="Remove all temporary files generated by the workflow. Use together "
        "with --dry-run to list files without actually deleting anything. Note "
        "that this will not recurse into subworkflows.",
    )
    group_utils.add_argument(
        "--bash-completion",
        action="store_true",
        help="Output code to register bash completion for snakemake. Put the "
        "following in your .bashrc (including the accents): "
        "`snakemake --bash-completion` or issue it in an open terminal "
        "session.",
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
        "Provenance-information based reports (e.g. --report and the "
        "--list_x_changes functions) will be empty or incomplete.",
    )
    group_utils.add_argument("--version", "-v", action="version", version=__version__)

    group_output = parser.add_argument_group("OUTPUT")
    group_output.add_argument(
        "--gui",
        nargs="?",
        const="8000",
        metavar="PORT",
        type=str,
        help="Serve an HTML based user interface to the given network and "
        "port e.g. 168.129.10.15:8000. By default Snakemake is only "
        "available in the local network (default port: 8000). To make "
        "Snakemake listen to all ip addresses add the special host address "
        "0.0.0.0 to the url (0.0.0.0:8000). This is important if Snakemake "
        "is used in a virtualised environment like Docker. If possible, a "
        "browser window is opened.",
    )
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
        "--stats",
        metavar="FILE",
        help="Write stats about Snakefile execution in JSON format to the given file.",
    )
    group_output.add_argument(
        "--nocolor", action="store_true", help="Do not use a colored output."
    )
    group_output.add_argument(
        "--quiet",
        "-q",
        nargs="*",
        choices=["progress", "rules", "all"],
        default=None,
        help="Do not output certain information. "
        "If used without arguments, do not output any progress or rule "
        "information. Defining 'all' results in no information being "
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
        "--latency-wait",
        "--output-wait",
        "-w",
        type=int,
        default=5,
        metavar="SECONDS",
        help="Wait given seconds if an output file of a job is not present after "
        "the job finished. This helps if your filesystem "
        "suffers from latency (default 5).",
    )
    group_behavior.add_argument(
        "--wait-for-files",
        nargs="*",
        metavar="FILE",
        help="Wait --latency-wait seconds for these "
        "files to be present before executing the workflow. "
        "This option is used internally to handle filesystem latency in cluster "
        "environments.",
    )
    group_behavior.add_argument(
        "--wait-for-files-file",
        metavar="FILE",
        help="Same behaviour as --wait-for-files, but file list is "
        "stored in file instead of being passed on the commandline. "
        "This is useful when the list of files is too long to be "
        "passed on the commandline.",
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
        "--keep-remote",
        action="store_true",
        help="Keep local copies of remote input files.",
    )
    group_behavior.add_argument(
        "--keep-target-files",
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
        "--target-jobs",
        nargs="+",
        help="Target particular jobs by RULE:WILDCARD1=VALUE,WILDCARD2=VALUE,... "
        "This is meant for internal use by Snakemake itself only.",
    )
    group_behavior.add_argument(
        "--local-groupid",
        default="local",
        help="Name for local groupid, meant for internal use only.",
    )
    group_behavior.add_argument(
        "--max-jobs-per-second",
        default=10,
        type=float,
        help="Maximal number of cluster/drmaa jobs per second, default is 10, "
        "fractions allowed.",
    )
    group_behavior.add_argument(
        "--max-status-checks-per-second",
        default=10,
        type=float,
        help="Maximal number of job status checks per second, default is 10, "
        "fractions allowed.",
    )
    group_behavior.add_argument(
        "-T",
        "--retries",
        "--restart-times",
        default=0,
        type=int,
        help="Number of times to restart failing jobs (defaults to 0).",
    )
    group_behavior.add_argument(
        "--attempt",
        default=1,
        type=int,
        help="Internal use only: define the initial value of the attempt "
        "parameter (default: 1).",
    )
    group_behavior.add_argument(
        "--wrapper-prefix",
        default="https://github.com/snakemake/snakemake-wrappers/raw/",
        help="Prefix for URL created from wrapper directive (default: "
        "https://github.com/snakemake/snakemake-wrappers/raw/). Set this to "
        "a different URL to use your fork or a local clone of the repository, "
        "e.g., use a git URL like 'git+file://path/to/your/local/clone@'.",
    )
    group_behavior.add_argument(
        "--default-remote-provider",
        choices=[
            "S3",
            "GS",
            "FTP",
            "SFTP",
            "S3Mocked",
            "gfal",
            "gridftp",
            "iRODS",
            "AzBlob",
            "XRootD",
        ],
        help="Specify default remote provider to be used for "
        "all input and output files that don't yet specify "
        "one.",
    )
    group_behavior.add_argument(
        "--default-remote-prefix",
        default="",
        help="Specify prefix for default remote provider. E.g. a bucket name.",
    )
    group_behavior.add_argument(
        "--no-shared-fs",
        action="store_true",
        help="Do not assume that jobs share a common file "
        "system. When this flag is activated, Snakemake will "
        "assume that the filesystem on a cluster node is not "
        "shared with other nodes. For example, this will lead "
        "to downloading remote files on each cluster node "
        "separately. Further, it won't take special measures "
        "to deal with filesystem latency issues. This option "
        "will in most cases only make sense in combination with "
        "--default-remote-provider. Further, when using --cluster "
        "you will have to also provide --cluster-status. "
        "Only activate this if you "
        "know what you are doing.",
    )
    group_behavior.add_argument(
        "--greediness",
        type=float,
        default=None,
        help="Set the greediness of scheduling. This value between 0 and 1 "
        "determines how careful jobs are selected for execution. The default "
        "value (1.0) provides the best speed and still acceptable scheduling "
        "quality.",
    )
    group_behavior.add_argument(
        "--no-hooks",
        action="store_true",
        help="Do not invoke onstart, onsuccess or onerror hooks after execution.",
    )
    group_behavior.add_argument(
        "--overwrite-shellcmd",
        help="Provide a shell command that shall be executed instead of those "
        "given in the workflow. "
        "This is for debugging purposes only.",
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
        "--mode",
        choices=[ExecMode.default, ExecMode.subprocess, ExecMode.remote],
        default=ExecMode.default,
        type=int,
        help="Set execution mode of Snakemake (internal use only).",
    )
    group_behavior.add_argument(
        "--show-failed-logs",
        action="store_true",
        help="Automatically display logs of failed jobs.",
    )
    group_behavior.add_argument(
        "--log-handler-script",
        metavar="FILE",
        default=None,
        help="Provide a custom script containing a function 'def log_handler(msg):'. "
        "Snakemake will call this function for every logging output (given as a dictionary msg)"
        "allowing to e.g. send notifications in the form of e.g. slack messages or emails.",
    )
    group_behavior.add_argument(
        "--log-service",
        default=None,
        choices=["none", "slack", "wms"],
        help="Set a specific messaging service for logging output."
        "Snakemake will notify the service on errors and completed execution."
        "Currently slack and workflow management system (wms) are supported.",
    )

    group_slurm = parser.add_argument_group("SLURM")
    slurm_mode_group = group_slurm.add_mutually_exclusive_group()

    slurm_mode_group.add_argument(
        "--slurm",
        action="store_true",
        help=(
            "Execute snakemake rules as SLURM batch jobs according"
            " to their 'resources' definition. SLURM resources as "
            " 'partition', 'ntasks', 'cpus', etc. need to be defined"
            " per rule within the 'resources' definition. Note, that"
            " memory can only be defined as 'mem_mb' or 'mem_mb_per_cpu'"
            " as analogous to the SLURM 'mem' and 'mem-per-cpu' flags"
            " to sbatch, respectively. Here, the unit is always 'MiB'."
            " In addition '--default_resources' should contain the"
            " SLURM account."
        ),
    ),
    slurm_mode_group.add_argument(
        "--slurm-jobstep",
        action="store_true",
        help=configargparse.SUPPRESS,  # this should be hidden and only be used
        # for snakemake to be working in jobscript-
        # mode
    )

    group_cluster = parser.add_argument_group("CLUSTER")

    # TODO extend below description to explain the wildcards that can be used
    cluster_mode_group = group_cluster.add_mutually_exclusive_group()
    cluster_mode_group.add_argument(
        "--cluster",
        metavar="CMD",
        help=(
            "Execute snakemake rules with the given submit command, "
            "e.g. qsub. Snakemake compiles jobs into scripts that are "
            "submitted to the cluster with the given command, once all input "
            "files for a particular job are present.\n"
            "The submit command can be decorated to make it aware of certain "
            "job properties (name, rulename, input, output, params, wildcards, log, threads "
            "and dependencies (see the argument below)), e.g.:\n"
            "$ snakemake --cluster 'qsub -pe threaded {threads}'."
        ),
    ),
    cluster_mode_group.add_argument(
        "--cluster-sync",
        metavar="CMD",
        help=(
            "cluster submission command will block, returning the remote exit"
            "status upon remote termination (for example, this should be used"
            "if the cluster command is 'qsub -sync y' (SGE)"
        ),
    ),
    cluster_mode_group.add_argument(
        "--drmaa",
        nargs="?",
        const="",
        metavar="ARGS",
        help="Execute snakemake on a cluster accessed via DRMAA, "
        "Snakemake compiles jobs into scripts that are "
        "submitted to the cluster with the given command, once all input "
        "files for a particular job are present. ARGS can be used to "
        "specify options of the underlying cluster system, "
        "thereby using the job properties name, rulename, input, output, params, wildcards, log, "
        "threads and dependencies, e.g.: "
        "--drmaa ' -pe threaded {threads}'. Note that ARGS must be given in quotes and "
        "with a leading whitespace.",
    )

    group_cluster.add_argument(
        "--immediate-submit",
        "--is",
        action="store_true",
        help="Immediately submit all jobs to the cluster instead of waiting "
        "for present input files. This will fail, unless you make "
        "the cluster aware of job dependencies, e.g. via:\n"
        "$ snakemake --cluster 'sbatch --dependency {dependencies}.\n"
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
        "The default script resides as 'jobscript.sh' in the "
        "installation directory.",
    )
    group_cluster.add_argument(
        "--jobname",
        "--jn",
        default="snakejob.{name}.{jobid}.sh",
        metavar="NAME",
        help="Provide a custom name for the jobscript that is submitted to the "
        'cluster (see --cluster). NAME is "snakejob.{name}.{jobid}.sh" '
        "per default. The wildcard {jobid} has to be present in the name.",
    )
    group_cluster.add_argument(
        "--cluster-status",
        help="Status command for cluster execution. This is only considered "
        "in combination with the --cluster flag. If provided, Snakemake will "
        "use the status command to determine if a job has finished successfully "
        "or failed. For this it is necessary that the submit command provided "
        "to --cluster returns the cluster job id. Then, the status command "
        "will be invoked with the job id. Snakemake expects it to return "
        "'success' if the job was successful, 'failed' if the job failed and "
        "'running' if the job still runs.",
    )
    group_cluster.add_argument(
        "--cluster-cancel",
        default=None,
        help="Specify a command that allows to stop currently running jobs. "
        "The command will be passed a single argument, the job id.",
    )
    group_cluster.add_argument(
        "--cluster-cancel-nargs",
        type=int,
        default=1000,
        help="Specify maximal number of job ids to pass to --cluster-cancel "
        "command, defaults to 1000.",
    )
    group_cluster.add_argument(
        "--cluster-sidecar",
        default=None,
        help="Optional command to start a sidecar process during cluster "
        "execution.  Only active when --cluster is given as well.",
    )
    group_cluster.add_argument(
        "--drmaa-log-dir",
        metavar="DIR",
        help="Specify a directory in which stdout and stderr files of DRMAA"
        " jobs will be written. The value may be given as a relative path,"
        " in which case Snakemake will use the current invocation directory"
        " as the origin. If given, this will override any given '-o' and/or"
        " '-e' native specification. If not given, all DRMAA stdout and"
        " stderr files are written to the current working directory.",
    )

    group_cloud = parser.add_argument_group("CLOUD")
    group_flux = parser.add_argument_group("FLUX")
    group_kubernetes = parser.add_argument_group("KUBERNETES")
    group_google_life_science = parser.add_argument_group("GOOGLE_LIFE_SCIENCE")
    group_kubernetes = parser.add_argument_group("KUBERNETES")
    group_tes = parser.add_argument_group("TES")
    group_tibanna = parser.add_argument_group("TIBANNA")

    group_kubernetes.add_argument(
        "--kubernetes",
        metavar="NAMESPACE",
        nargs="?",
        const="default",
        help="Execute workflow in a kubernetes cluster (in the cloud). "
        "NAMESPACE is the namespace you want to use for your job (if nothing "
        "specified: 'default'). "
        "Usually, this requires --default-remote-provider and "
        "--default-remote-prefix to be set to a S3 or GS bucket where your . "
        "data shall be stored. It is further advisable to activate conda "
        "integration via --use-conda.",
    )
    group_kubernetes.add_argument(
        "--container-image",
        metavar="IMAGE",
        help="Docker image to use, e.g., when submitting jobs to kubernetes "
        "Defaults to 'https://hub.docker.com/r/snakemake/snakemake', tagged with "
        "the same version as the currently running Snakemake instance. "
        "Note that overwriting this value is up to your responsibility. "
        "Any used image has to contain a working snakemake installation "
        "that is compatible with (or ideally the same as) the currently "
        "running version.",
    )
    group_kubernetes.add_argument(
        "--k8s-cpu-scalar",
        metavar="FLOAT",
        default=0.95,
        type=float,
        help="K8s reserves some proportion of available CPUs for its own use. "
        "So, where an underlying node may have 8 CPUs, only e.g. 7600 milliCPUs "
        "are allocatable to k8s pods (i.e. snakemake jobs). As 8 > 7.6, k8s can't "
        "find a node with enough CPU resource to run such jobs. This argument acts "
        "as a global scalar on each job's CPU request, so that e.g. a job whose "
        "rule definition asks for 8 CPUs will request 7600m CPUs from k8s, "
        "allowing it to utilise one entire node. N.B: the job itself would still "
        "see the original value, i.e. as the value substituted in {threads}.",
    )

    group_tibanna.add_argument(
        "--tibanna",
        action="store_true",
        help="Execute workflow on AWS cloud using Tibanna. This requires "
        "--default-remote-prefix to be set to S3 bucket name and prefix"
        " (e.g. 'bucketname/subdirectory') where input is already stored"
        " and output will be sent to. Using --tibanna implies --default-resources"
        " is set as default. Optionally, use --precommand to"
        " specify any preparation command to run before snakemake command"
        " on the cloud (inside snakemake container on Tibanna VM)."
        " Also, --use-conda, --use-singularity, --config, --configfile are"
        " supported and will be carried over.",
    )
    group_tibanna.add_argument(
        "--tibanna-sfn",
        help="Name of Tibanna Unicorn step function (e.g. tibanna_unicorn_monty)."
        "This works as serverless scheduler/resource allocator and must be "
        "deployed first using tibanna cli. (e.g. tibanna deploy_unicorn --usergroup="
        "monty --buckets=bucketname)",
    )
    group_tibanna.add_argument(
        "--precommand",
        help="Any command to execute before snakemake command on AWS cloud "
        "such as wget, git clone, unzip, etc. This is used with --tibanna."
        "Do not include input/output download/upload commands - file transfer"
        " between S3 bucket and the run environment (container) is automatically"
        " handled by Tibanna.",
    )
    group_tibanna.add_argument(
        "--tibanna-config",
        nargs="+",
        help="Additional tibanna config e.g. --tibanna-config spot_instance=true subnet="
        "<subnet_id> security group=<security_group_id>",
    )
    group_google_life_science.add_argument(
        "--google-lifesciences",
        action="store_true",
        help="Execute workflow on Google Cloud cloud using the Google Life. "
        " Science API. This requires default application credentials (json) "
        " to be created and export to the environment to use Google Cloud "
        " Storage, Compute Engine, and Life Sciences. The credential file "
        " should be exported as GOOGLE_APPLICATION_CREDENTIALS for snakemake "
        " to discover. Also, --use-conda, --use-singularity, --config, "
        "--configfile are supported and will be carried over.",
    )
    group_google_life_science.add_argument(
        "--google-lifesciences-regions",
        nargs="+",
        default=["us-east1", "us-west1", "us-central1"],
        help="Specify one or more valid instance regions (defaults to US)",
    )
    group_google_life_science.add_argument(
        "--google-lifesciences-location",
        help="The Life Sciences API service used to schedule the jobs. "
        " E.g., us-centra1 (Iowa) and europe-west2 (London) "
        " Watch the terminal output to see all options found to be available. "
        " If not specified, defaults to the first found with a matching prefix "
        " from regions specified with --google-lifesciences-regions.",
    )
    group_google_life_science.add_argument(
        "--google-lifesciences-keep-cache",
        action="store_true",
        help="Cache workflows in your Google Cloud Storage Bucket specified "
        "by --default-remote-prefix/{source}/{cache}. Each workflow working "
        "directory is compressed to a .tar.gz, named by the hash of the "
        "contents, and kept in Google Cloud Storage. By default, the caches "
        "are deleted at the shutdown step of the workflow.",
    )

    group_azure_batch = parser.add_argument_group("AZURE_BATCH")

    group_azure_batch.add_argument(
        "--az-batch",
        action="store_true",
        help="Execute workflow on azure batch",
    )

    group_azure_batch.add_argument(
        "--az-batch-enable-autoscale",
        action="store_true",
        help="Enable autoscaling of the azure batch pool nodes, this option will set the initial dedicated node count to zero, and requires five minutes to resize the cluster, so is only recommended for longer running jobs.",
    )

    group_azure_batch.add_argument(
        "--az-batch-account-url",
        nargs="?",
        help="Azure batch account url, requires AZ_BATCH_ACCOUNT_KEY environment variable to be set.",
    )

    group_flux.add_argument(
        "--flux",
        action="store_true",
        help="Execute your workflow on a flux cluster. "
        "Flux can work with both a shared network filesystem (like NFS) or without. "
        "If you don't have a shared filesystem, additionally specify --no-shared-fs.",
    )

    group_tes.add_argument(
        "--tes",
        metavar="URL",
        help="Send workflow tasks to GA4GH TES server specified by url.",
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
        help="Specify a directory in which the 'conda' and 'conda-archive' "
        "directories are created. These are used to store conda environments "
        "and their archives, respectively. If not supplied, the value is set "
        "to the '.snakemake' directory relative to the invocation directory. "
        "If supplied, the `--use-conda` flag must also be set. The value may "
        "be given as a relative path, which will be extrapolated to the "
        "invocation directory, or as an absolute path. The value can also be "
        "provided via the environment variable $SNAKEMAKE_CONDA_PREFIX.",
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
        nargs="?",
        help="Cleanup conda packages after creating environments. "
        "In case of 'tarballs' mode, will clean up all downloaded package tarballs. "
        "In case of 'cache' mode, will additionally clean up unused package caches. "
        "If mode is omitted, will default to only cleaning up the tarballs.",
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
        default="mamba",
        choices=["conda", "mamba"],
        help="Choose the conda frontend for installing environments. "
        "Mamba is much faster and highly recommended.",
    )

    group_singularity = parser.add_argument_group("SINGULARITY")

    group_singularity.add_argument(
        "--use-singularity",
        action="store_true",
        help="If defined in the rule, run job within a singularity container. "
        "If this flag is not set, the singularity directive is ignored.",
    )
    group_singularity.add_argument(
        "--singularity-prefix",
        metavar="DIR",
        help="Specify a directory in which singularity images will be stored."
        "If not supplied, the value is set "
        "to the '.snakemake' directory relative to the invocation directory. "
        "If supplied, the `--use-singularity` flag must also be set. The value "
        "may be given as a relative path, which will be extrapolated to the "
        "invocation directory, or as an absolute path.",
    )
    group_singularity.add_argument(
        "--singularity-args",
        default="",
        metavar="ARGS",
        help="Pass additional args to singularity.",
    )
    group_singularity.add_argument(
        "--cleanup-containers",
        action="store_true",
        help="Remove unused (singularity) containers",
    )

    group_env_modules = parser.add_argument_group("ENVIRONMENT MODULES")

    group_env_modules.add_argument(
        "--use-envmodules",
        action="store_true",
        help="If defined in the rule, run job within the given environment "
        "modules, loaded in the given order. This can be combined with "
        "--use-conda and --use-singularity, which will then be only used as a "
        "fallback for rules which don't define environment modules.",
    )

    # Add namespaced arguments to parser for each plugin
    ExecutorPluginRegistry().register_cli_args(parser)
    return parser


def generate_parser_metadata(parser, args):
    """Given a populated parser, generate the original command along with
    metadata that can be handed to a logger to use as needed.
    """
    command = "snakemake %s" % " ".join(
        parser._source_to_settings["command_line"][""][1]
    )
    workdir = os.getcwd()
    metadata = args.__dict__
    metadata.update({"command": command})
    return metadata


def main(argv=None):
    """Main entry point."""

    if sys.version_info < MIN_PY_VERSION:
        print(
            f"Snakemake requires at least Python {MIN_PY_VERSION}.",
            file=sys.stderr,
        )
        exit(1)

    parser = get_argument_parser()
    args = parser.parse_args(argv)

    snakefile = args.snakefile
    if snakefile is None:
        for p in SNAKEFILE_CHOICES:
            if os.path.exists(p):
                snakefile = p
                break
        if snakefile is None:
            print(
                "Error: no Snakefile found, tried {}.".format(
                    ", ".join(SNAKEFILE_CHOICES)
                ),
                file=sys.stderr,
            )
            sys.exit(1)

    # Custom argument parsing based on chosen executor
    # We also only validate an executor plugin when it's selected
    executor_args = None
    if args.executor:
        plugin = ExecutorPluginRegistry().plugins[args.executor]

        # This is the dataclass prepared by the executor
        executor_args = plugin.get_executor_settings(args)

        # Hold a handle to the plugin class
        executor_args._executor = plugin

    workflow_profile = None
    if args.workflow_profile != "none":
        if args.workflow_profile:
            workflow_profile = args.workflow_profile
        else:
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

    if (args.profile or workflow_profile) and args.mode == ExecMode.default:
        # Reparse args while inferring config file from profile.
        # But only do this if the user has invoked Snakemake (ExecMode.default)
        profiles = []
        if args.profile:
            profiles.append(args.profile)
        if workflow_profile:
            profiles.append(workflow_profile)

        print(
            f"Using profile{'s' if len(profiles) > 1 else ''} "
            f"{' and '.join(map(str, profiles))} for setting default command line arguments.",
            file=sys.stderr,
        )

        parser = get_argument_parser(profiles=profiles)
        args = parser.parse_args(argv)

        def adjust_path(f):
            if os.path.exists(f) or os.path.isabs(f):
                return f
            else:
                return get_profile_file(args.profile, f, return_default=True)

        # update file paths to be relative to the profile
        # (if they do not exist relative to CWD)
        if args.jobscript:
            args.jobscript = adjust_path(args.jobscript)
        if args.cluster:
            args.cluster = adjust_path(args.cluster)
        if args.cluster_sync:
            args.cluster_sync = adjust_path(args.cluster_sync)
        for key in "cluster_status", "cluster_cancel", "cluster_sidecar":
            if getattr(args, key):
                setattr(args, key, adjust_path(getattr(args, key)))
        if args.report_stylesheet:
            args.report_stylesheet = adjust_path(args.report_stylesheet)

    if args.quiet is not None and len(args.quiet) == 0:
        # default case, set quiet to progress and rule
        args.quiet = ["progress", "rules"]

    if args.bash_completion:
        cmd = b"complete -o bashdefault -C snakemake-bash-completion snakemake"
        sys.stdout.buffer.write(cmd)
        sys.exit(0)

    if args.batch is not None and args.forceall:
        print(
            "--batch may not be combined with --forceall, because recomputed upstream "
            "jobs in subsequent batches may render already obtained results outdated."
        )

    try:
        resources = parse_resources(args.resources)
        config = parse_config(args)

        if args.default_resources is not None:
            default_resources = DefaultResources(args.default_resources)
        else:
            default_resources = None

        batch = parse_batch(args)
        overwrite_threads = parse_set_threads(args)
        overwrite_resources = parse_set_resources(args)
        overwrite_resource_scopes = parse_set_resource_scope(args)

        overwrite_scatter = parse_set_scatter(args)

        overwrite_groups = parse_groups(args)
        group_components = parse_group_components(args)
    except ValueError as e:
        print(e, file=sys.stderr)
        print("", file=sys.stderr)
        sys.exit(1)

    non_local_exec = (
        args.cluster
        or args.slurm
        or args.slurm_jobstep
        or args.cluster_sync
        or args.tibanna
        or args.kubernetes
        or args.tes
        or args.az_batch
        or args.google_lifesciences
        or args.drmaa
        or args.flux
    )
    no_exec = (
        args.print_compilation
        or args.list_code_changes
        or args.list_conda_envs
        or args.list_input_changes
        or args.list_params_changes
        or args.list
        or args.list_target_rules
        or args.list_untracked
        or args.list_version_changes
        or args.export_cwl
        or args.generate_unit_tests
        or args.dag
        or args.d3dag
        or args.filegraph
        or args.rulegraph
        or args.summary
        or args.detailed_summary
        or args.lint
        or args.containerize
        or args.report
        or args.gui
        or args.archive
        or args.unlock
        or args.cleanup_metadata
    )

    try:
        cores, jobs = parse_cores_jobs(
            args.cores, args.jobs, no_exec, non_local_exec, args.dryrun
        )
        args.cores = cores
        args.jobs = jobs
    except CliException as err:
        print(err.msg, sys.stderr)
        sys.exit(1)

    if args.drmaa_log_dir is not None:
        if not os.path.isabs(args.drmaa_log_dir):
            args.drmaa_log_dir = os.path.abspath(os.path.expanduser(args.drmaa_log_dir))

    if args.runtime_profile:
        import yappi

        yappi.start()

    if args.immediate_submit and not args.notemp:
        print(
            "Error: --immediate-submit has to be combined with --notemp, "
            "because temp file handling is not supported in this mode.",
            file=sys.stderr,
        )
        sys.exit(1)

    if (args.conda_prefix or args.conda_create_envs_only) and not args.use_conda:
        if args.conda_prefix and os.environ.get("SNAKEMAKE_CONDA_PREFIX", False):
            print(
                "Warning: The enviorment variable SNAKEMAKE_CONDA_PREFIX is set"
                "but --use-conda is not."
                "Snakemake will ignore SNAKEMAKE_CONDA_PREFIX"
                "and conda enviorments will not be used or created.",
                file=sys.stderr,
            )
            args.conda_prefix = None
        else:
            print(
                "Error: --use-conda must be set if --conda-prefix or "
                "--create-envs-only is set.",
                file=sys.stderr,
            )
            sys.exit(1)

    if args.singularity_prefix and not args.use_singularity:
        print(
            "Error: --use_singularity must be set if --singularity-prefix is set.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.kubernetes and (
        not args.default_remote_provider or not args.default_remote_prefix
    ):
        print(
            "Error: --kubernetes must be combined with "
            "--default-remote-provider and --default-remote-prefix, see "
            "https://snakemake.readthedocs.io/en/stable/executing/cloud.html"
            "#executing-a-snakemake-workflow-via-kubernetes",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.tibanna:
        if not args.default_remote_prefix:
            print(
                "Error: --tibanna must be combined with --default-remote-prefix "
                "to provide bucket name and subdirectory (prefix) "
                "(e.g. 'bucketname/projectname'",
                file=sys.stderr,
            )
            sys.exit(1)
        args.default_remote_prefix = args.default_remote_prefix.rstrip("/")
        if not args.tibanna_sfn:
            args.tibanna_sfn = os.environ.get("TIBANNA_DEFAULT_STEP_FUNCTION_NAME", "")
            if not args.tibanna_sfn:
                print(
                    "Error: to use --tibanna, either --tibanna-sfn or environment variable "
                    "TIBANNA_DEFAULT_STEP_FUNCTION_NAME must be set and exported "
                    "to provide name of the tibanna unicorn step function "
                    "(e.g. 'tibanna_unicorn_monty'). The step function must be deployed first "
                    "using tibanna cli (e.g. tibanna deploy_unicorn --usergroup=monty "
                    "--buckets=bucketname)",
                    file=sys.stderr,
                )
                sys.exit(1)

    if args.az_batch:
        if not args.default_remote_provider or not args.default_remote_prefix:
            print(
                "Error: --az-batch must be combined with "
                "--default-remote-provider AzBlob and --default-remote-prefix to "
                "provide a blob container name\n",
                file=sys.stderr,
            )
            sys.exit(1)
        elif args.az_batch_account_url is None:
            print(
                "Error: --az-batch-account-url must be set when --az-batch is used\n",
                file=sys.stderr,
            )
            sys.exit(1)
        elif not url_can_parse(args.az_batch_account_url):
            print(
                "Error: invalide azure batch account url, please use format: https://{account_name}.{location}.batch.azure.com."
            )
            sys.exit(1)
        elif os.getenv("AZ_BATCH_ACCOUNT_KEY") is None:
            print(
                "Error: environment variable AZ_BATCH_ACCOUNT_KEY must be set when --az-batch is used\n",
                file=sys.stderr,
            )
            sys.exit(1)

    if args.google_lifesciences:
        if not os.environ.get("GOOGLE_APPLICATION_CREDENTIALS"):
            print(
                "Error: GOOGLE_APPLICATION_CREDENTIALS environment variable must "
                "be available for --google-lifesciences",
                file=sys.stderr,
            )
            sys.exit(1)

        if not args.default_remote_prefix:
            print(
                "Error: --google-life-sciences must be combined with "
                " --default-remote-prefix to provide bucket name and "
                "subdirectory (prefix) (e.g. 'bucketname/projectname'",
                file=sys.stderr,
            )
            sys.exit(1)

    if args.delete_all_output and args.delete_temp_output:
        print(
            "Error: --delete-all-output and --delete-temp-output are mutually exclusive.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.gui is not None:
        try:
            import snakemake.gui as gui
        except ImportError:
            print(
                "Error: GUI needs Flask to be installed. Install "
                "with easy_install or contact your administrator.",
                file=sys.stderr,
            )
            sys.exit(1)

        _logging.getLogger("werkzeug").setLevel(_logging.ERROR)

        _snakemake = partial(snakemake, os.path.abspath(snakefile))
        gui.register(_snakemake, args)

        if ":" in args.gui:
            host, port = args.gui.split(":")
        else:
            port = args.gui
            host = "127.0.0.1"

        url = f"http://{host}:{port}"
        print(f"Listening on {url}.", file=sys.stderr)

        def open_browser():
            try:
                webbrowser.open(url)
            except:
                pass

        print("Open this address in your browser to access the GUI.", file=sys.stderr)
        threading.Timer(0.5, open_browser).start()
        success = True

        try:
            gui.app.run(debug=False, threaded=True, port=int(port), host=host)

        except (KeyboardInterrupt, SystemExit):
            # silently close
            pass
    else:
        log_handler = []
        if args.log_handler_script is not None:
            if not os.path.exists(args.log_handler_script):
                print(
                    "Error: no log handler script found, {}.".format(
                        args.log_handler_script
                    ),
                    file=sys.stderr,
                )
                sys.exit(1)
            log_script = SourceFileLoader("log", args.log_handler_script).load_module()
            try:
                log_handler.append(log_script.log_handler)
            except:
                print(
                    'Error: Invalid log handler script, {}. Expect python function "log_handler(msg)".'.format(
                        args.log_handler_script
                    ),
                    file=sys.stderr,
                )
                sys.exit(1)

        if args.log_service == "slack":
            slack_logger = logging.SlackLogger()
            log_handler.append(slack_logger.log_handler)

        elif args.wms_monitor or args.log_service == "wms":
            # Generate additional metadata for server
            metadata = generate_parser_metadata(parser, args)
            wms_logger = logging.WMSLogger(
                args.wms_monitor, args.wms_monitor_arg, metadata=metadata
            )
            log_handler.append(wms_logger.log_handler)

        if args.draft_notebook:
            from snakemake import notebook

            args.target = [args.draft_notebook]
            args.edit_notebook = notebook.EditMode(draft_only=True)
        elif args.edit_notebook:
            from snakemake import notebook

            args.target = [args.edit_notebook]
            args.force = True
            args.edit_notebook = notebook.EditMode(args.notebook_listen)

        aggregated_wait_for_files = args.wait_for_files
        if args.wait_for_files_file is not None:
            wait_for_files([args.wait_for_files_file], latency_wait=args.latency_wait)

            with open(args.wait_for_files_file) as fd:
                extra_wait_files = [line.strip() for line in fd.readlines()]

            if aggregated_wait_for_files is None:
                aggregated_wait_for_files = extra_wait_files
            else:
                aggregated_wait_for_files.extend(extra_wait_files)

        success = snakemake(
            snakefile,
            batch=batch,
            cache=args.cache,
            report=args.report,
            report_stylesheet=args.report_stylesheet,
            lint=args.lint,
            containerize=args.containerize,
            generate_unit_tests=args.generate_unit_tests,
            listrules=args.list,
            list_target_rules=args.list_target_rules,
            cores=args.cores,
            local_cores=args.local_cores,
            nodes=args.jobs,
            resources=resources,
            overwrite_threads=overwrite_threads,
            max_threads=args.max_threads,
            overwrite_scatter=overwrite_scatter,
            default_resources=default_resources,
            overwrite_resources=overwrite_resources,
            overwrite_resource_scopes=overwrite_resource_scopes,
            config=config,
            configfiles=args.configfile,
            config_args=args.config,
            workdir=args.directory,
            targets=args.target,
            target_jobs=parse_target_jobs_cli_args(args),
            dryrun=args.dryrun,
            printshellcmds=args.printshellcmds,
            debug_dag=args.debug_dag,
            printdag=args.dag,
            printrulegraph=args.rulegraph,
            printfilegraph=args.filegraph,
            printd3dag=args.d3dag,
            touch=args.touch,
            forcetargets=args.force,
            forceall=args.forceall,
            forcerun=args.forcerun,
            prioritytargets=args.prioritize,
            until=args.until,
            omit_from=args.omit_from,
            stats=args.stats,
            nocolor=args.nocolor,
            quiet=args.quiet,
            keepgoing=args.keep_going,
            slurm=args.slurm,
            slurm_jobstep=args.slurm_jobstep,
            rerun_triggers=args.rerun_triggers,
            cluster=args.cluster,
            cluster_sync=args.cluster_sync,
            drmaa=args.drmaa,
            drmaa_log_dir=args.drmaa_log_dir,
            kubernetes=args.kubernetes,
            container_image=args.container_image,
            k8s_cpu_scalar=args.k8s_cpu_scalar,
            flux=args.flux,
            tibanna=args.tibanna,
            tibanna_sfn=args.tibanna_sfn,
            az_batch=args.az_batch,
            az_batch_enable_autoscale=args.az_batch_enable_autoscale,
            az_batch_account_url=args.az_batch_account_url,
            google_lifesciences=args.google_lifesciences,
            google_lifesciences_regions=args.google_lifesciences_regions,
            google_lifesciences_location=args.google_lifesciences_location,
            google_lifesciences_cache=args.google_lifesciences_keep_cache,
            tes=args.tes,
            precommand=args.precommand,
            preemption_default=args.preemption_default,
            preemptible_rules=args.preemptible_rules,
            tibanna_config=args.tibanna_config,
            jobname=args.jobname,
            immediate_submit=args.immediate_submit,
            standalone=True,
            ignore_ambiguity=args.allow_ambiguity,
            lock=not args.nolock,
            unlock=args.unlock,
            cleanup_metadata=args.cleanup_metadata,
            conda_cleanup_envs=args.conda_cleanup_envs,
            cleanup_containers=args.cleanup_containers,
            cleanup_shadow=args.cleanup_shadow,
            force_incomplete=args.rerun_incomplete,
            ignore_incomplete=args.ignore_incomplete,
            list_version_changes=args.list_version_changes,
            list_code_changes=args.list_code_changes,
            list_input_changes=args.list_input_changes,
            list_params_changes=args.list_params_changes,
            list_untracked=args.list_untracked,
            summary=args.summary,
            detailed_summary=args.detailed_summary,
            archive=args.archive,
            delete_all_output=args.delete_all_output,
            delete_temp_output=args.delete_temp_output,
            print_compilation=args.print_compilation,
            verbose=args.verbose,
            debug=args.debug,
            jobscript=args.jobscript,
            notemp=args.notemp,
            all_temp=args.all_temp,
            keep_remote_local=args.keep_remote,
            greediness=args.greediness,
            no_hooks=args.no_hooks,
            overwrite_shellcmd=args.overwrite_shellcmd,
            latency_wait=args.latency_wait,
            wait_for_files=aggregated_wait_for_files,
            keep_target_files=args.keep_target_files,
            allowed_rules=args.allowed_rules,
            max_jobs_per_second=args.max_jobs_per_second,
            max_status_checks_per_second=args.max_status_checks_per_second,
            restart_times=args.retries,
            attempt=args.attempt,
            force_use_threads=args.force_use_threads,
            use_conda=args.use_conda,
            conda_frontend=args.conda_frontend,
            conda_prefix=args.conda_prefix,
            conda_cleanup_pkgs=args.conda_cleanup_pkgs,
            list_conda_envs=args.list_conda_envs,
            use_singularity=args.use_singularity,
            use_env_modules=args.use_envmodules,
            singularity_prefix=args.singularity_prefix,
            shadow_prefix=args.shadow_prefix,
            singularity_args=args.singularity_args,
            scheduler=args.scheduler,
            scheduler_ilp_solver=args.scheduler_ilp_solver,
            conda_create_envs_only=args.conda_create_envs_only,
            mode=args.mode,
            wrapper_prefix=args.wrapper_prefix,
            default_remote_provider=args.default_remote_provider,
            default_remote_prefix=args.default_remote_prefix,
            assume_shared_fs=not args.no_shared_fs,
            cluster_status=args.cluster_status,
            cluster_cancel=args.cluster_cancel,
            cluster_cancel_nargs=args.cluster_cancel_nargs,
            cluster_sidecar=args.cluster_sidecar,
            export_cwl=args.export_cwl,
            show_failed_logs=args.show_failed_logs,
            keep_incomplete=args.keep_incomplete,
            keep_metadata=not args.drop_metadata,
            edit_notebook=args.edit_notebook,
            envvars=args.envvars,
            overwrite_groups=overwrite_groups,
            group_components=group_components,
            max_inventory_wait_time=args.max_inventory_time,
            log_handler=log_handler,
            execute_subworkflows=not args.no_subworkflows,
            conda_not_block_search_path_envvars=args.conda_not_block_search_path_envvars,
            scheduler_solver_path=args.scheduler_solver_path,
            conda_base_path=args.conda_base_path,
            local_groupid=args.local_groupid,
            executor_args=executor_args,
            cleanup_scripts=not args.skip_script_cleanup,
        )

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

    sys.exit(0 if success else 1)


def bash_completion(snakefile="Snakefile"):
    """Entry point for bash completion."""
    if not len(sys.argv) >= 2:
        print(
            "Calculate bash completion for snakemake. This tool shall not be invoked by hand."
        )
        sys.exit(1)

    def print_candidates(candidates):
        if candidates:
            candidates = sorted(set(candidates))
            ## Use bytes for avoiding '^M' under Windows.
            sys.stdout.buffer.write(b"\n".join(s.encode() for s in candidates))

    prefix = sys.argv[2]

    if prefix.startswith("-"):
        print_candidates(
            action.option_strings[0]
            for action in get_argument_parser()._actions
            if action.option_strings and action.option_strings[0].startswith(prefix)
        )
    else:
        candidates = []
        files = glob.glob(f"{prefix}*")
        if files:
            candidates.extend(files)
        if os.path.exists(snakefile):
            workflow = Workflow(snakefile=snakefile)
            workflow.include(snakefile)

            candidates.extend(
                [file for file in workflow.concrete_files if file.startswith(prefix)]
                + [rule.name for rule in workflow.rules if rule.name.startswith(prefix)]
            )
        if len(candidates) > 0:
            print_candidates(candidates)
    sys.exit(0)