import os

from snakemake.shell import shell
from snakemake.logging import logger
from snakemake.exceptions import WorkflowError


def subworkflow(
    path,
    basedir,
    input,
    output,
    params,
    wildcards,
    threads,
    resources,
    log,
    config,
    rulename,
    conda_env,
    container_img,
    singularity_args,
    env_modules,
    bench_record,
    jobid,
    bench_iteration,
    cleanup_scripts,
    shadow_dir,
):
    # TODO: use `subsnakemake` instead?

    # extract (and normalize) parameters
    basedir = os.path.normpath(basedir)
    sub_workdir = os.path.normpath(
        params.get("workdir", f"subworkflow_{rulename}")
    )
    sub_configfile = params.get("configfile", "")

    snakemake_params = params.get("snakemake_params", "")

    # the sub-workflow needs to run in a different directory than the
    # main-workflow. This is necessary due to Snakemake's locking mechanism
    if basedir == sub_workdir:
        raise WorkflowError(
            "Sub-workflow workdir must be different from main-workflow workdir"
        )

    # all output files must use the sub-workflow's workdir as their prefix
    # otherwise file dependencies cannot be easily handled
    output_normalized = []
    for fname in output:
        if not fname.startswith(sub_workdir):
            raise WorkflowError(
                "All output files must begin with "
                f"(normalized) sub-workflow workdir ({sub_workdir})"
            )

        fname_rel = os.path.relpath(fname, start=sub_workdir)
        output_normalized.append(fname_rel)

    # prepare parameters
    configfile_arg = ""
    if len(sub_configfile) > 0:
        configfile_arg = f"--configfile {sub_configfile}"

    # execute sub-workflow
    logger.info("Executing sub-workflow")

    shell(
        "snakemake"
        "  --snakefile {path}"
        "  --directory {sub_workdir}"
        "  --cores {threads}"
        "  {configfile_arg}"
        "  -F {snakemake_params}"
        "  {output_normalized}"
    )

    logger.info("Finished execution of sub-workflow")
