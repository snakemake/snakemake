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
    raise WorkflowError(
        "Rule-like subworkflows not yet implemented."
    )
