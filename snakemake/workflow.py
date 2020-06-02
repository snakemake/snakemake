__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2019, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import re
import os
import sys
import signal
import json
import urllib
from collections import OrderedDict
from itertools import filterfalse, chain
from functools import partial
from operator import attrgetter
import copy
import subprocess
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import pathname2url, url2pathname

from snakemake.logging import logger, format_resources, format_resource_names
from snakemake.rules import Rule, Ruleorder, RuleProxy
from snakemake.exceptions import (
    RuleException,
    CreateRuleException,
    UnknownRuleException,
    NoRulesException,
    print_exception,
    WorkflowError,
)
from snakemake.shell import shell
from snakemake.dag import DAG
from snakemake.scheduler import JobScheduler
from snakemake.parser import parse
import snakemake.io
from snakemake.io import (
    protected,
    temp,
    temporary,
    ancient,
    directory,
    expand,
    dynamic,
    glob_wildcards,
    flag,
    not_iterable,
    touch,
    unpack,
    local,
    pipe,
    repeat,
    report,
    multiext,
    IOFile,
)
from snakemake.persistence import Persistence
from snakemake.utils import update_config
from snakemake.script import script
from snakemake.notebook import notebook
from snakemake.wrapper import wrapper
from snakemake.cwl import cwl
import snakemake.wrapper
from snakemake.common import Mode, bytesto, ON_WINDOWS
from snakemake.utils import simplify_path
from snakemake.checkpoints import Checkpoint, Checkpoints
from snakemake.resources import DefaultResources
from snakemake.caching.local import OutputFileCache as LocalOutputFileCache
from snakemake.caching.remote import OutputFileCache as RemoteOutputFileCache


class Workflow:
    def __init__(
        self,
        snakefile=None,
        jobscript=None,
        overwrite_shellcmd=None,
        overwrite_config=dict(),
        overwrite_workdir=None,
        overwrite_configfiles=None,
        overwrite_clusterconfig=dict(),
        overwrite_threads=dict(),
        config_args=None,
        debug=False,
        verbose=False,
        use_conda=False,
        conda_frontend=None,
        conda_prefix=None,
        use_singularity=False,
        use_env_modules=False,
        singularity_prefix=None,
        singularity_args="",
        shadow_prefix=None,
        mode=Mode.default,
        wrapper_prefix=None,
        printshellcmds=False,
        restart_times=None,
        attempt=1,
        default_remote_provider=None,
        default_remote_prefix="",
        run_local=True,
        default_resources=None,
        cache=None,
        nodes=1,
        cores=1,
        resources=None,
        conda_cleanup_pkgs=None,
        edit_notebook=False,
        envvars=None,
    ):
        """
        Create the controller.
        """

        self.global_resources = dict() if resources is None else resources
        self.global_resources["_cores"] = cores
        self.global_resources["_nodes"] = nodes

        self._rules = OrderedDict()
        self.first_rule = None
        self._workdir = None
        self.overwrite_workdir = overwrite_workdir
        self.workdir_init = os.path.abspath(os.curdir)
        self._ruleorder = Ruleorder()
        self._localrules = set()
        self.linemaps = dict()
        self.rule_count = 0
        self.basedir = os.path.dirname(snakefile)
        self.snakefile = os.path.abspath(snakefile)
        self.included = []
        self.included_stack = []
        self.jobscript = jobscript
        self.persistence = None
        self.globals = globals()
        self._subworkflows = dict()
        self.overwrite_shellcmd = overwrite_shellcmd
        self.overwrite_config = overwrite_config
        self.overwrite_configfiles = overwrite_configfiles
        self.overwrite_clusterconfig = overwrite_clusterconfig
        self.overwrite_threads = overwrite_threads
        self.config_args = config_args
        self.immediate_submit = None
        self._onsuccess = lambda log: None
        self._onerror = lambda log: None
        self._onstart = lambda log: None
        self._wildcard_constraints = dict()
        self.debug = debug
        self.verbose = verbose
        self._rulecount = 0
        self.use_conda = use_conda
        self.conda_frontend = conda_frontend
        self.conda_prefix = conda_prefix
        self.use_singularity = use_singularity
        self.use_env_modules = use_env_modules
        self.singularity_prefix = singularity_prefix
        self.singularity_args = singularity_args
        self.shadow_prefix = shadow_prefix
        self.global_container_img = None
        self.mode = mode
        self.wrapper_prefix = wrapper_prefix
        self.printshellcmds = printshellcmds
        self.restart_times = restart_times
        self.attempt = attempt
        self.default_remote_provider = default_remote_provider
        self.default_remote_prefix = default_remote_prefix
        self.configfiles = []
        self.run_local = run_local
        self.report_text = None
        self.conda_cleanup_pkgs = conda_cleanup_pkgs
        self.edit_notebook = edit_notebook
        # environment variables to pass to jobs
        # These are defined via the "envvars:" syntax in the Snakefile itself
        self.envvars = set()

        self.enable_cache = False
        if cache is not None:
            self.enable_cache = True
            self.cache_rules = set(cache)
            if self.default_remote_provider is not None:
                self.output_file_cache = RemoteOutputFileCache(
                    self.default_remote_provider
                )
            else:
                self.output_file_cache = LocalOutputFileCache()
        else:
            self.output_file_cache = None
            self.cache_rules = set()

        if default_resources is not None:
            self.default_resources = default_resources
        else:
            # only _cores and _nodes
            self.default_resources = DefaultResources()

        self.iocache = snakemake.io.IOCache()

        global config
        config = copy.deepcopy(self.overwrite_config)

        global cluster_config
        cluster_config = copy.deepcopy(self.overwrite_clusterconfig)

        global rules
        rules = Rules()
        global checkpoints
        checkpoints = Checkpoints()

        if envvars is not None:
            self.register_envvars(*envvars)

    def lint(self, json=False):
        from snakemake.linting.rules import RuleLinter
        from snakemake.linting.snakefiles import SnakefileLinter

        json_snakefile_lints, snakefile_linted = SnakefileLinter(
            self, self.included
        ).lint(json=json)
        json_rule_lints, rules_linted = RuleLinter(self, self.rules).lint(json=json)

        linted = snakefile_linted or rules_linted

        if json:
            import json

            print(
                json.dumps(
                    {"snakefiles": json_snakefile_lints, "rules": json_rule_lints},
                    indent=2,
                )
            )
        else:
            if not linted:
                logger.info("Congratulations, your workflow is in a good condition!")
        return linted

    def is_cached_rule(self, rule: Rule):
        return rule.name in self.cache_rules

    def get_sources(self):
        files = set()

        def local_path(f):
            if ON_WINDOWS:
                try:
                    f = pathname2url(f)
                except OSError:
                    pass  # f isn't changed if it wasn't a path

            url = urlparse(f)
            if url.scheme == "file" or url.scheme == "":
                return url2pathname(url.path)
            return None

        def norm_rule_relpath(f, rule):
            if not os.path.isabs(f):
                f = os.path.join(rule.basedir, f)
            return os.path.relpath(f)

        # get registered sources
        for f in self.included:
            f = local_path(f)
            if f:
                try:
                    f = os.path.relpath(f)
                except ValueError:
                    if ON_WINDOWS:
                        pass  # relpath doesn't work on win if files are on different drive
                    else:
                        raise
                files.add(f)
        for rule in self.rules:
            script_path = rule.script or rule.notebook
            if script_path:
                script_path = norm_rule_relpath(script_path, rule)
                files.add(script_path)
                script_dir = os.path.dirname(script_path)
                files.update(
                    os.path.join(dirpath, f)
                    for dirpath, _, files in os.walk(script_dir)
                    for f in files
                )
            if rule.conda_env:
                f = local_path(rule.conda_env)
                if f:
                    # url points to a local env file
                    env_path = norm_rule_relpath(f, rule)
                    files.add(env_path)

        for f in self.configfiles:
            files.add(f)

        # get git-managed files
        # TODO allow a manifest file as alternative
        try:
            out = subprocess.check_output(
                ["git", "ls-files", "--recurse-submodules", "."], stderr=subprocess.PIPE
            )
            for f in out.decode().split("\n"):
                if f:
                    files.add(os.path.relpath(f))
        except subprocess.CalledProcessError as e:
            if "fatal: not a git repository" in e.stderr.decode().lower():
                logger.warning(
                    "Unable to retrieve additional files from git. "
                    "This is not a git repository."
                )
            else:
                raise WorkflowError(
                    "Error executing git:\n{}".format(e.stderr.decode())
                )

        return files

    def check_source_sizes(self, filename, warning_size_gb=0.2):
        """A helper function to check the filesize, and return the file
           to the calling function Additionally, given that we encourage these 
           packages to be small, we set a warning at 200MB (0.2GB).
        """
        gb = bytesto(os.stat(filename).st_size, "g")
        if gb > warning_size_gb:
            logger.warning(
                "File {} (size {} GB) is greater than the {} GB suggested size "
                "Consider uploading larger files to storage first.".format(
                    filename, gb, warning_size_gb
                )
            )
        return filename

    @property
    def subworkflows(self):
        return self._subworkflows.values()

    @property
    def rules(self):
        return self._rules.values()

    @property
    def cores(self):
        return self.global_resources["_cores"]

    @property
    def nodes(self):
        return self.global_resources["_nodes"]

    @property
    def concrete_files(self):
        return (
            file
            for rule in self.rules
            for file in chain(rule.input, rule.output)
            if not callable(file) and not file.contains_wildcard()
        )

    def check(self):
        for clause in self._ruleorder:
            for rulename in clause:
                if not self.is_rule(rulename):
                    raise UnknownRuleException(
                        rulename, prefix="Error in ruleorder definition."
                    )

    def add_rule(self, name=None, lineno=None, snakefile=None, checkpoint=False):
        """
        Add a rule.
        """
        if name is None:
            name = str(len(self._rules) + 1)
        if self.is_rule(name):
            raise CreateRuleException(
                "The name {} is already used by another rule".format(name)
            )
        rule = Rule(name, self, lineno=lineno, snakefile=snakefile)
        self._rules[rule.name] = rule
        self.rule_count += 1
        if not self.first_rule:
            self.first_rule = rule.name
        return name

    def is_rule(self, name):
        """
        Return True if name is the name of a rule.

        Arguments
        name -- a name
        """
        return name in self._rules

    def get_rule(self, name):
        """
        Get rule by name.

        Arguments
        name -- the name of the rule
        """
        if not self._rules:
            raise NoRulesException()
        if not name in self._rules:
            raise UnknownRuleException(name)
        return self._rules[name]

    def list_rules(self, only_targets=False):
        rules = self.rules
        if only_targets:
            rules = filterfalse(Rule.has_wildcards, rules)
        for rule in rules:
            logger.rule_info(name=rule.name, docstring=rule.docstring)

    def list_resources(self):
        for resource in set(
            resource for rule in self.rules for resource in rule.resources
        ):
            if resource not in "_cores _nodes".split():
                logger.info(resource)

    def is_local(self, rule):
        return rule.group is None and (rule.name in self._localrules or rule.norun)

    def check_localrules(self):
        undefined = self._localrules - set(rule.name for rule in self.rules)
        if undefined:
            logger.warning(
                "localrules directive specifies rules that are not "
                "present in the Snakefile:\n{}\n".format(
                    "\n".join(map("\t{}".format, undefined))
                )
            )

    def inputfile(self, path):
        """Mark file as being an input file of the workflow.

        This also means that eventual --default-remote-provider/prefix settings
        will be applied to this file. The file is returned as _IOFile object,
        such that it can e.g. be transparently opened with _IOFile.open().
        """
        if isinstance(path, Path):
            path = str(path)
        if self.default_remote_provider is not None:
            path = self.apply_default_remote(path)
        return IOFile(path)

    def apply_default_remote(self, path):
        """Apply the defined default remote provider to the given path and return the updated _IOFile.
        Asserts that default remote provider is defined.
        """
        assert (
            self.default_remote_provider is not None
        ), "No default remote provider is defined, calling this anyway is a bug"
        path = "{}/{}".format(self.default_remote_prefix, path)
        path = os.path.normpath(path)
        return self.default_remote_provider.remote(path)

    def execute(
        self,
        targets=None,
        dryrun=False,
        touch=False,
        local_cores=1,
        forcetargets=False,
        forceall=False,
        forcerun=None,
        until=[],
        omit_from=[],
        prioritytargets=None,
        quiet=False,
        keepgoing=False,
        printshellcmds=False,
        printreason=False,
        printdag=False,
        cluster=None,
        cluster_sync=None,
        jobname=None,
        immediate_submit=False,
        ignore_ambiguity=False,
        printrulegraph=False,
        printfilegraph=False,
        printd3dag=False,
        drmaa=None,
        drmaa_log_dir=None,
        kubernetes=None,
        tibanna=None,
        tibanna_sfn=None,
        google_lifesciences=None,
        google_lifesciences_regions=None,
        google_lifesciences_location=None,
        google_lifesciences_cache=False,
        precommand="",
        tibanna_config=False,
        container_image=None,
        stats=None,
        force_incomplete=False,
        ignore_incomplete=False,
        list_version_changes=False,
        list_code_changes=False,
        list_input_changes=False,
        list_params_changes=False,
        list_untracked=False,
        list_conda_envs=False,
        summary=False,
        archive=None,
        delete_all_output=False,
        delete_temp_output=False,
        detailed_summary=False,
        latency_wait=3,
        wait_for_files=None,
        nolock=False,
        unlock=False,
        notemp=False,
        nodeps=False,
        cleanup_metadata=None,
        conda_cleanup_envs=False,
        cleanup_shadow=False,
        cleanup_scripts=True,
        subsnakemake=None,
        updated_files=None,
        keep_target_files=False,
        keep_shadow=False,
        keep_remote_local=False,
        allowed_rules=None,
        max_jobs_per_second=None,
        max_status_checks_per_second=None,
        greediness=1.0,
        no_hooks=False,
        force_use_threads=False,
        conda_create_envs_only=False,
        assume_shared_fs=True,
        cluster_status=None,
        report=None,
        report_stylesheet=None,
        export_cwl=False,
        batch=None,
        keepincomplete=False,
    ):

        self.check_localrules()
        self.immediate_submit = immediate_submit
        self.cleanup_scripts = cleanup_scripts

        def rules(items):
            return map(self._rules.__getitem__, filter(self.is_rule, items))

        if keep_target_files:

            def files(items):
                return filterfalse(self.is_rule, items)

        else:

            def files(items):
                relpath = (
                    lambda f: f
                    if os.path.isabs(f) or f.startswith("root://")
                    else os.path.relpath(f)
                )
                return map(relpath, filterfalse(self.is_rule, items))

        if not targets:
            targets = [self.first_rule] if self.first_rule is not None else list()

        if prioritytargets is None:
            prioritytargets = list()
        if forcerun is None:
            forcerun = list()
        if until is None:
            until = list()
        if omit_from is None:
            omit_from = list()

        priorityrules = set(rules(prioritytargets))
        priorityfiles = set(files(prioritytargets))
        forcerules = set(rules(forcerun))
        forcefiles = set(files(forcerun))
        untilrules = set(rules(until))
        untilfiles = set(files(until))
        omitrules = set(rules(omit_from))
        omitfiles = set(files(omit_from))
        targetrules = set(
            chain(
                rules(targets),
                filterfalse(Rule.has_wildcards, priorityrules),
                filterfalse(Rule.has_wildcards, forcerules),
                filterfalse(Rule.has_wildcards, untilrules),
            )
        )
        targetfiles = set(chain(files(targets), priorityfiles, forcefiles, untilfiles))
        if forcetargets:
            forcefiles.update(targetfiles)
            forcerules.update(targetrules)

        rules = self.rules
        if allowed_rules:
            rules = [rule for rule in rules if rule.name in set(allowed_rules)]

        if wait_for_files is not None:
            try:
                snakemake.io.wait_for_files(wait_for_files, latency_wait=latency_wait)
            except IOError as e:
                logger.error(str(e))
                return False

        dag = DAG(
            self,
            rules,
            dryrun=dryrun,
            targetfiles=targetfiles,
            targetrules=targetrules,
            # when cleaning up conda, we should enforce all possible jobs
            # since their envs shall not be deleted
            forceall=forceall or conda_cleanup_envs,
            forcefiles=forcefiles,
            forcerules=forcerules,
            priorityfiles=priorityfiles,
            priorityrules=priorityrules,
            untilfiles=untilfiles,
            untilrules=untilrules,
            omitfiles=omitfiles,
            omitrules=omitrules,
            ignore_ambiguity=ignore_ambiguity,
            force_incomplete=force_incomplete,
            ignore_incomplete=ignore_incomplete
            or printdag
            or printrulegraph
            or printfilegraph,
            notemp=notemp,
            keep_remote_local=keep_remote_local,
            batch=batch,
        )

        self.persistence = Persistence(
            nolock=nolock,
            dag=dag,
            conda_prefix=self.conda_prefix,
            singularity_prefix=self.singularity_prefix,
            shadow_prefix=self.shadow_prefix,
            warn_only=dryrun
            or printrulegraph
            or printfilegraph
            or printdag
            or summary
            or archive
            or list_version_changes
            or list_code_changes
            or list_input_changes
            or list_params_changes
            or list_untracked
            or delete_all_output
            or delete_temp_output,
        )

        if cleanup_metadata:
            for f in cleanup_metadata:
                self.persistence.cleanup_metadata(f)
            return True

        logger.info("Building DAG of jobs...")
        dag.init()
        dag.update_checkpoint_dependencies()
        # check incomplete has to run BEFORE any call to postprocess
        dag.check_incomplete()
        dag.check_dynamic()

        if unlock:
            try:
                self.persistence.cleanup_locks()
                logger.info("Unlocking working directory.")
                return True
            except IOError:
                logger.error(
                    "Error: Unlocking the directory {} failed. Maybe "
                    "you don't have the permissions?"
                )
                return False
        try:
            self.persistence.lock()
        except IOError:
            logger.error(
                "Error: Directory cannot be locked. Please make "
                "sure that no other Snakemake process is trying to create "
                "the same files in the following directory:\n{}\n"
                "If you are sure that no other "
                "instances of snakemake are running on this directory, "
                "the remaining lock was likely caused by a kill signal or "
                "a power loss. It can be removed with "
                "the --unlock argument.".format(os.getcwd())
            )
            return False

        if cleanup_shadow:
            self.persistence.cleanup_shadow()
            return True

        if (
            self.subworkflows
            and not printdag
            and not printrulegraph
            and not printfilegraph
        ):
            # backup globals
            globals_backup = dict(self.globals)
            # execute subworkflows
            for subworkflow in self.subworkflows:
                subworkflow_targets = subworkflow.targets(dag)
                logger.debug(
                    "Files requested from subworkflow:\n    {}".format(
                        "\n    ".join(subworkflow_targets)
                    )
                )
                updated = list()
                if subworkflow_targets:
                    logger.info("Executing subworkflow {}.".format(subworkflow.name))
                    if not subsnakemake(
                        subworkflow.snakefile,
                        workdir=subworkflow.workdir,
                        targets=subworkflow_targets,
                        configfiles=[subworkflow.configfile]
                        if subworkflow.configfile
                        else None,
                        updated_files=updated,
                    ):
                        return False
                    dag.updated_subworkflow_files.update(
                        subworkflow.target(f) for f in updated
                    )
                else:
                    logger.info(
                        "Subworkflow {}: Nothing to be done.".format(subworkflow.name)
                    )
            if self.subworkflows:
                logger.info("Executing main workflow.")
            # rescue globals
            self.globals.update(globals_backup)

        dag.postprocess()
        # deactivate IOCache such that from now on we always get updated
        # size, existence and mtime information
        # ATTENTION: this may never be removed without really good reason.
        # Otherwise weird things may happen.
        self.iocache.deactivate()
        # clear and deactivate persistence cache, from now on we want to see updates
        self.persistence.deactivate_cache()

        if nodeps:
            missing_input = [
                f
                for job in dag.targetjobs
                for f in job.input
                if dag.needrun(job) and not os.path.exists(f)
            ]
            if missing_input:
                logger.error(
                    "Dependency resolution disabled (--nodeps) "
                    "but missing input "
                    "files detected. If this happens on a cluster, please make sure "
                    "that you handle the dependencies yourself or turn off "
                    "--immediate-submit. Missing input files:\n{}".format(
                        "\n".join(missing_input)
                    )
                )
                return False

        updated_files.extend(f for job in dag.needrun_jobs for f in job.output)

        if export_cwl:
            from snakemake.cwl import dag_to_cwl
            import json

            with open(export_cwl, "w") as cwl:
                json.dump(dag_to_cwl(dag), cwl, indent=4)
            return True
        elif report:
            from snakemake.report import auto_report

            auto_report(dag, report, stylesheet=report_stylesheet)
            return True
        elif printd3dag:
            dag.d3dag()
            return True
        elif printdag:
            print(dag)
            return True
        elif printrulegraph:
            print(dag.rule_dot())
            return True
        elif printfilegraph:
            print(dag.filegraph_dot())
            return True
        elif summary:
            print("\n".join(dag.summary(detailed=False)))
            return True
        elif detailed_summary:
            print("\n".join(dag.summary(detailed=True)))
            return True
        elif archive:
            dag.archive(archive)
            return True
        elif delete_all_output:
            dag.clean(only_temp=False, dryrun=dryrun)
            return True
        elif delete_temp_output:
            dag.clean(only_temp=True, dryrun=dryrun)
            return True
        elif list_version_changes:
            items = list(chain(*map(self.persistence.version_changed, dag.jobs)))
            if items:
                print(*items, sep="\n")
            return True
        elif list_code_changes:
            items = list(chain(*map(self.persistence.code_changed, dag.jobs)))
            for j in dag.jobs:
                items.extend(list(j.outputs_older_than_script_or_notebook()))
            if items:
                print(*items, sep="\n")
            return True
        elif list_input_changes:
            items = list(chain(*map(self.persistence.input_changed, dag.jobs)))
            if items:
                print(*items, sep="\n")
            return True
        elif list_params_changes:
            items = list(chain(*map(self.persistence.params_changed, dag.jobs)))
            if items:
                print(*items, sep="\n")
            return True
        elif list_untracked:
            dag.list_untracked()
            return True

        if self.use_singularity:
            if assume_shared_fs:
                dag.pull_container_imgs(
                    dryrun=dryrun or list_conda_envs, quiet=list_conda_envs
                )
        if self.use_conda:
            if assume_shared_fs:
                dag.create_conda_envs(
                    dryrun=dryrun or list_conda_envs or conda_cleanup_envs,
                    quiet=list_conda_envs,
                )
            if conda_create_envs_only:
                return True

        if list_conda_envs:
            print("environment", "container", "location", sep="\t")
            for env in set(job.conda_env for job in dag.jobs):
                if env:
                    print(
                        simplify_path(env.file),
                        env.container_img_url or "",
                        simplify_path(env.path),
                        sep="\t",
                    )
            return True

        if conda_cleanup_envs:
            self.persistence.conda_cleanup_envs()
            return True

        scheduler = JobScheduler(
            self,
            dag,
            self.cores,
            local_cores=local_cores,
            dryrun=dryrun,
            touch=touch,
            cluster=cluster,
            cluster_status=cluster_status,
            cluster_config=cluster_config,
            cluster_sync=cluster_sync,
            jobname=jobname,
            max_jobs_per_second=max_jobs_per_second,
            max_status_checks_per_second=max_status_checks_per_second,
            quiet=quiet,
            keepgoing=keepgoing,
            drmaa=drmaa,
            drmaa_log_dir=drmaa_log_dir,
            kubernetes=kubernetes,
            tibanna=tibanna,
            tibanna_sfn=tibanna_sfn,
            google_lifesciences=google_lifesciences,
            google_lifesciences_regions=google_lifesciences_regions,
            google_lifesciences_location=google_lifesciences_location,
            google_lifesciences_cache=google_lifesciences_cache,
            precommand=precommand,
            tibanna_config=tibanna_config,
            container_image=container_image,
            printreason=printreason,
            printshellcmds=printshellcmds,
            latency_wait=latency_wait,
            greediness=greediness,
            force_use_threads=force_use_threads,
            assume_shared_fs=assume_shared_fs,
            keepincomplete=keepincomplete,
        )

        if not dryrun:
            if len(dag):
                shell_exec = shell.get_executable()
                if shell_exec is not None:
                    logger.info("Using shell: {}".format(shell_exec))
                if cluster or cluster_sync or drmaa:
                    logger.resources_info(
                        "Provided cluster nodes: {}".format(self.nodes)
                    )
                elif kubernetes or tibanna:
                    logger.resources_info("Provided cloud nodes: {}".format(self.nodes))
                else:
                    if self.cores is not None:
                        warning = (
                            ""
                            if self.cores > 1
                            else " (use --cores to define parallelism)"
                        )
                        logger.resources_info(
                            "Provided cores: {}{}".format(self.cores, warning)
                        )
                        logger.resources_info(
                            "Rules claiming more threads " "will be scaled down."
                        )

                provided_resources = format_resources(self.global_resources)
                if provided_resources:
                    logger.resources_info("Provided resources: " + provided_resources)

                if self.run_local and any(rule.group for rule in self.rules):
                    logger.info("Group jobs: inactive (local execution)")

                if not self.use_conda and any(rule.conda_env for rule in self.rules):
                    logger.info("Conda environments: ignored")

                if not self.use_singularity and any(
                    rule.container_img for rule in self.rules
                ):
                    logger.info("Singularity containers: ignored")

                logger.run_info("\n".join(dag.stats()))
            else:
                logger.info("Nothing to be done.")
        else:
            # the dryrun case
            if len(dag):
                logger.run_info("\n".join(dag.stats()))
            else:
                logger.info("Nothing to be done.")
                return True
            if quiet:
                # in case of dryrun and quiet, just print above info and exit
                return True

        if not dryrun and not no_hooks:
            self._onstart(logger.get_logfile())

        success = scheduler.schedule()

        if not immediate_submit:
            dag.cleanup_workdir()

        if success:
            if dryrun:
                if len(dag):
                    logger.run_info("\n".join(dag.stats()))
                logger.info(
                    "This was a dry-run (flag -n). The order of jobs "
                    "does not reflect the order of execution."
                )
                logger.remove_logfile()
            else:
                if stats:
                    scheduler.stats.to_json(stats)
                logger.logfile_hint()
            if not dryrun and not no_hooks:
                self._onsuccess(logger.get_logfile())
            return True
        else:
            if not dryrun and not no_hooks:
                self._onerror(logger.get_logfile())
            logger.logfile_hint()
            return False

    @property
    def current_basedir(self):
        """Basedir of currently parsed Snakefile."""
        assert self.included_stack
        return os.path.abspath(os.path.dirname(self.included_stack[-1]))

    def register_envvars(self, *envvars):
        """
        Register environment variables that shall be passed to jobs.
        If used multiple times, union is taken.
        """
        undefined = set(var for var in envvars if var not in os.environ)
        if undefined:
            raise WorkflowError(
                "The following environment variables are requested by the workflow but undefined. "
                "Please make sure that they are correctly defined before running Snakemake:\n"
                "{}".format("\n".join(undefined))
            )
        self.envvars.update(envvars)

    def include(
        self,
        snakefile,
        overwrite_first_rule=False,
        print_compilation=False,
        overwrite_shellcmd=None,
    ):
        """
        Include a snakefile.
        """
        if isinstance(snakefile, Path):
            snakefile = str(snakefile)

        # check if snakefile is a path to the filesystem
        if not urllib.parse.urlparse(snakefile).scheme:
            if not os.path.isabs(snakefile) and self.included_stack:
                snakefile = os.path.join(self.current_basedir, snakefile)
            # Could still be an url if relative import was used
            if not urllib.parse.urlparse(snakefile).scheme:
                snakefile = os.path.abspath(snakefile)
        # else it could be an url.
        # at least we don't want to modify the path for clarity.

        if snakefile in self.included:
            logger.info("Multiple include of {} ignored".format(snakefile))
            return
        self.included.append(snakefile)
        self.included_stack.append(snakefile)

        global workflow

        workflow = self

        first_rule = self.first_rule
        code, linemap, rulecount = parse(
            snakefile,
            overwrite_shellcmd=self.overwrite_shellcmd,
            rulecount=self._rulecount,
        )
        self._rulecount = rulecount

        if print_compilation:
            print(code)

        # insert the current directory into sys.path
        # this allows to import modules from the workflow directory
        sys.path.insert(0, os.path.dirname(snakefile))

        self.linemaps[snakefile] = linemap
        exec(compile(code, snakefile, "exec"), self.globals)
        if not overwrite_first_rule:
            self.first_rule = first_rule
        self.included_stack.pop()

    def onstart(self, func):
        """Register onstart function."""
        self._onstart = func

    def onsuccess(self, func):
        """Register onsuccess function."""
        self._onsuccess = func

    def onerror(self, func):
        """Register onerror function."""
        self._onerror = func

    def global_wildcard_constraints(self, **content):
        """Register global wildcard constraints."""
        self._wildcard_constraints.update(content)
        # update all rules so far
        for rule in self.rules:
            rule.update_wildcard_constraints()

    def workdir(self, workdir):
        """Register workdir."""
        if self.overwrite_workdir is None:
            os.makedirs(workdir, exist_ok=True)
            self._workdir = workdir
            os.chdir(workdir)

    def configfile(self, fp):
        """ Update the global config with data from the given file. """
        global config
        self.configfiles.append(fp)
        c = snakemake.io.load_configfile(fp)
        update_config(config, c)
        update_config(config, self.overwrite_config)

    def report(self, path):
        """ Define a global report description in .rst format."""
        self.report_text = os.path.join(self.current_basedir, path)

    @property
    def config(self):
        global config
        return config

    def ruleorder(self, *rulenames):
        self._ruleorder.add(*rulenames)

    def subworkflow(self, name, snakefile=None, workdir=None, configfile=None):
        # Take absolute path of config file, because it is relative to current
        # workdir, which could be changed for the subworkflow.
        if configfile:
            configfile = os.path.abspath(configfile)
        sw = Subworkflow(self, name, snakefile, workdir, configfile)
        self._subworkflows[name] = sw
        self.globals[name] = sw.target

    def localrules(self, *rulenames):
        self._localrules.update(rulenames)

    def rule(self, name=None, lineno=None, snakefile=None, checkpoint=False):
        name = self.add_rule(name, lineno, snakefile, checkpoint)
        rule = self.get_rule(name)
        rule.is_checkpoint = checkpoint

        def decorate(ruleinfo):
            if ruleinfo.wildcard_constraints:
                rule.set_wildcard_constraints(
                    *ruleinfo.wildcard_constraints[0],
                    **ruleinfo.wildcard_constraints[1]
                )
            if ruleinfo.input:
                rule.set_input(*ruleinfo.input[0], **ruleinfo.input[1])
            if ruleinfo.output:
                rule.set_output(*ruleinfo.output[0], **ruleinfo.output[1])
            if ruleinfo.params:
                rule.set_params(*ruleinfo.params[0], **ruleinfo.params[1])
            # handle default resources
            if self.default_resources is not None:
                rule.resources = copy.deepcopy(self.default_resources.parsed)
            if ruleinfo.threads is not None:
                if (
                    not isinstance(ruleinfo.threads, int)
                    and not isinstance(ruleinfo.threads, float)
                    and not callable(ruleinfo.threads)
                ):
                    raise RuleException(
                        "Threads value has to be an integer, float, or a callable.",
                        rule=rule,
                    )
                if name in self.overwrite_threads:
                    rule.resources["_cores"] = self.overwrite_threads[name]
                else:
                    if isinstance(ruleinfo.threads, float):
                        ruleinfo.threads = int(ruleinfo.threads)
                    rule.resources["_cores"] = ruleinfo.threads
            if ruleinfo.shadow_depth:
                if ruleinfo.shadow_depth not in (True, "shallow", "full", "minimal"):
                    raise RuleException(
                        "Shadow must either be 'minimal', 'shallow', 'full', "
                        "or True (equivalent to 'full')",
                        rule=rule,
                    )
                if ruleinfo.shadow_depth is True:
                    rule.shadow_depth = "full"
                    logger.warning(
                        "Shadow is set to True in rule {} (equivalent to 'full'). It's encouraged to use the more explicit options 'minimal|shallow|full' instead.".format(
                            rule
                        )
                    )
                else:
                    rule.shadow_depth = ruleinfo.shadow_depth
            if ruleinfo.resources:
                args, resources = ruleinfo.resources
                if args:
                    raise RuleException("Resources have to be named.")
                if not all(
                    map(
                        lambda r: isinstance(r, int)
                        or isinstance(r, str)
                        or callable(r),
                        resources.values(),
                    )
                ):
                    raise RuleException(
                        "Resources values have to be integers, strings, or callables (functions)",
                        rule=rule,
                    )
                rule.resources.update(resources)
            if ruleinfo.priority:
                if not isinstance(ruleinfo.priority, int) and not isinstance(
                    ruleinfo.priority, float
                ):
                    raise RuleException(
                        "Priority values have to be numeric.", rule=rule
                    )
                rule.priority = ruleinfo.priority
            if ruleinfo.version:
                rule.version = ruleinfo.version
            if ruleinfo.log:
                rule.set_log(*ruleinfo.log[0], **ruleinfo.log[1])
            if ruleinfo.message:
                rule.message = ruleinfo.message
            if ruleinfo.benchmark:
                rule.benchmark = ruleinfo.benchmark
            if not self.run_local and ruleinfo.group is not None:
                rule.group = ruleinfo.group
            if ruleinfo.wrapper:
                rule.conda_env = snakemake.wrapper.get_conda_env(
                    ruleinfo.wrapper, prefix=self.wrapper_prefix
                )
                # TODO retrieve suitable singularity image

            if ruleinfo.env_modules:
                # If using environment modules and they are defined for the rule,
                # ignore conda and singularity directive below.
                # The reason is that this is likely intended in order to use
                # a software stack specifically compiled for a particular
                # HPC cluster.
                invalid_rule = not (
                    ruleinfo.script
                    or ruleinfo.wrapper
                    or ruleinfo.shellcmd
                    or ruleinfo.notebook
                )
                if invalid_rule:
                    raise RuleException(
                        "envmodules directive is only allowed with "
                        "shell, script, notebook, or wrapper directives (not with run)",
                        rule=rule,
                    )
                from snakemake.deployment.env_modules import EnvModules

                rule.env_modules = EnvModules(*ruleinfo.env_modules)

            if ruleinfo.conda_env:
                if not (
                    ruleinfo.script
                    or ruleinfo.wrapper
                    or ruleinfo.shellcmd
                    or ruleinfo.notebook
                ):
                    raise RuleException(
                        "Conda environments are only allowed "
                        "with shell, script, notebook, or wrapper directives "
                        "(not with run).",
                        rule=rule,
                    )
                if not (
                    urllib.parse.urlparse(ruleinfo.conda_env).scheme
                    or os.path.isabs(ruleinfo.conda_env)
                ):
                    ruleinfo.conda_env = os.path.join(
                        self.current_basedir, ruleinfo.conda_env
                    )
                rule.conda_env = ruleinfo.conda_env

            invalid_rule = not (
                ruleinfo.script
                or ruleinfo.wrapper
                or ruleinfo.shellcmd
                or ruleinfo.notebook
            )
            if ruleinfo.container_img:
                if invalid_rule:
                    raise RuleException(
                        "Singularity directive is only allowed "
                        "with shell, script, notebook or wrapper directives "
                        "(not with run).",
                        rule=rule,
                    )
                rule.container_img = ruleinfo.container_img
            elif self.global_container_img:
                if not invalid_rule:
                    # skip rules with run directive
                    rule.container_img = self.global_container_img

            rule.norun = ruleinfo.norun
            rule.docstring = ruleinfo.docstring
            rule.run_func = ruleinfo.func
            rule.shellcmd = ruleinfo.shellcmd
            rule.script = ruleinfo.script
            rule.notebook = ruleinfo.notebook
            rule.wrapper = ruleinfo.wrapper
            rule.cwl = ruleinfo.cwl
            rule.restart_times = self.restart_times
            rule.basedir = self.current_basedir

            if ruleinfo.cache is True:
                if not self.enable_cache:
                    logger.warning(
                        "Workflow defines that rule {} is eligible for caching between workflows "
                        "(use the --cache argument to enable this).".format(rule.name)
                    )
                else:
                    self.cache_rules.add(rule.name)
            elif not (ruleinfo.cache is False):
                raise WorkflowError(
                    "Invalid argument for 'cache:' directive. Only true allowed. "
                    "To deactivate caching, remove directive.",
                    rule=rule,
                )

            ruleinfo.func.__name__ = "__{}".format(rule.name)
            self.globals[ruleinfo.func.__name__] = ruleinfo.func
            setattr(rules, rule.name, RuleProxy(rule))
            if checkpoint:
                checkpoints.register(rule)
            return ruleinfo.func

        return decorate

    def docstring(self, string):
        def decorate(ruleinfo):
            ruleinfo.docstring = string
            return ruleinfo

        return decorate

    def input(self, *paths, **kwpaths):
        def decorate(ruleinfo):
            ruleinfo.input = (paths, kwpaths)
            return ruleinfo

        return decorate

    def output(self, *paths, **kwpaths):
        def decorate(ruleinfo):
            ruleinfo.output = (paths, kwpaths)
            return ruleinfo

        return decorate

    def params(self, *params, **kwparams):
        def decorate(ruleinfo):
            ruleinfo.params = (params, kwparams)
            return ruleinfo

        return decorate

    def wildcard_constraints(self, *wildcard_constraints, **kwwildcard_constraints):
        def decorate(ruleinfo):
            ruleinfo.wildcard_constraints = (
                wildcard_constraints,
                kwwildcard_constraints,
            )
            return ruleinfo

        return decorate

    def cache_rule(self, cache):
        def decorate(ruleinfo):
            ruleinfo.cache = cache
            return ruleinfo

        return decorate

    def message(self, message):
        def decorate(ruleinfo):
            ruleinfo.message = message
            return ruleinfo

        return decorate

    def benchmark(self, benchmark):
        def decorate(ruleinfo):
            ruleinfo.benchmark = benchmark
            return ruleinfo

        return decorate

    def conda(self, conda_env):
        def decorate(ruleinfo):
            ruleinfo.conda_env = conda_env
            return ruleinfo

        return decorate

    def container(self, container_img):
        def decorate(ruleinfo):
            ruleinfo.container_img = container_img
            return ruleinfo

        return decorate

    def envmodules(self, *env_modules):
        def decorate(ruleinfo):
            ruleinfo.env_modules = env_modules
            return ruleinfo

        return decorate

    def global_container(self, container_img):
        self.global_container_img = container_img

    def threads(self, threads):
        def decorate(ruleinfo):
            ruleinfo.threads = threads
            return ruleinfo

        return decorate

    def shadow(self, shadow_depth):
        def decorate(ruleinfo):
            ruleinfo.shadow_depth = shadow_depth
            return ruleinfo

        return decorate

    def resources(self, *args, **resources):
        def decorate(ruleinfo):
            ruleinfo.resources = (args, resources)
            return ruleinfo

        return decorate

    def priority(self, priority):
        def decorate(ruleinfo):
            ruleinfo.priority = priority
            return ruleinfo

        return decorate

    def version(self, version):
        def decorate(ruleinfo):
            ruleinfo.version = version
            return ruleinfo

        return decorate

    def group(self, group):
        def decorate(ruleinfo):
            ruleinfo.group = group
            return ruleinfo

        return decorate

    def log(self, *logs, **kwlogs):
        def decorate(ruleinfo):
            ruleinfo.log = (logs, kwlogs)
            return ruleinfo

        return decorate

    def shellcmd(self, cmd):
        def decorate(ruleinfo):
            ruleinfo.shellcmd = cmd
            return ruleinfo

        return decorate

    def script(self, script):
        def decorate(ruleinfo):
            ruleinfo.script = script
            return ruleinfo

        return decorate

    def notebook(self, notebook):
        def decorate(ruleinfo):
            ruleinfo.notebook = notebook
            return ruleinfo

        return decorate

    def wrapper(self, wrapper):
        def decorate(ruleinfo):
            ruleinfo.wrapper = wrapper
            return ruleinfo

        return decorate

    def cwl(self, cwl):
        def decorate(ruleinfo):
            ruleinfo.cwl = cwl
            return ruleinfo

        return decorate

    def norun(self):
        def decorate(ruleinfo):
            ruleinfo.norun = True
            return ruleinfo

        return decorate

    def run(self, func):
        return RuleInfo(func)

    @staticmethod
    def _empty_decorator(f):
        return f


class RuleInfo:
    def __init__(self, func):
        self.func = func
        self.shellcmd = None
        self.norun = False
        self.input = None
        self.output = None
        self.params = None
        self.message = None
        self.benchmark = None
        self.conda_env = None
        self.container_img = None
        self.env_modules = None
        self.wildcard_constraints = None
        self.threads = None
        self.shadow_depth = None
        self.resources = None
        self.priority = None
        self.version = None
        self.log = None
        self.docstring = None
        self.group = None
        self.script = None
        self.notebook = None
        self.wrapper = None
        self.cwl = None
        self.cache = False


class Subworkflow:
    def __init__(self, workflow, name, snakefile, workdir, configfile):
        self.workflow = workflow
        self.name = name
        self._snakefile = snakefile
        self._workdir = workdir
        self.configfile = configfile

    @property
    def snakefile(self):
        if self._snakefile is None:
            return os.path.abspath(os.path.join(self.workdir, "Snakefile"))
        if not os.path.isabs(self._snakefile):
            return os.path.abspath(os.path.join(self.workflow.basedir, self._snakefile))
        return self._snakefile

    @property
    def workdir(self):
        workdir = "." if self._workdir is None else self._workdir
        if not os.path.isabs(workdir):
            return os.path.abspath(os.path.join(self.workflow.basedir, workdir))
        return workdir

    def target(self, paths):
        if not_iterable(paths):
            path = paths
            path = (
                path
                if os.path.isabs(path) or path.startswith("root://")
                else os.path.join(self.workdir, path)
            )
            return flag(path, "subworkflow", self)
        return [self.target(path) for path in paths]

    def targets(self, dag):
        def relpath(f):
            if f.startswith(self.workdir):
                return os.path.relpath(f, start=self.workdir)
            # do not adjust absolute targets outside of workdir
            return f

        return [
            relpath(f)
            for job in dag.jobs
            for f in job.subworkflow_input
            if job.subworkflow_input[f] is self
        ]


class Rules:
    """ A namespace for rules so that they can be accessed via dot notation. """

    pass


def srcdir(path):
    """Return the absolute path, relative to the source directory of the current Snakefile."""
    if not workflow.included_stack:
        return None
    return os.path.join(os.path.dirname(workflow.included_stack[-1]), path)
