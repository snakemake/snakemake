import shutil
from itertools import chain

from snakemake.exceptions import WorkflowError
from snakemake.shell import shell


# class Nextflow:
#     def __init__(
#         self, 
#         input,
#         params,
#         threads,
#         resources,
#         log,
#         use_singularity,
#         use_conda,
#         bench_record,
#         is_pipeline,
#     ):
#         self.input = input
#         self.params = params
#         self.threads = threads
#         self.resources = resources
#         self.log = log
#         self.use_singularity = use_singularity
#         self.use_conda = use_conda
#         self.bench_record = bench_record
#         self.is_pipeline = is_pipeline

#     def run(self):


def nextflow(
    repo,
    input,
    params,
    threads,
    resources,
    log,
    use_singularity,
    use_conda,
    bench_record,
    is_pipeline,
):
    if shutil.which("nextflow") is None:
        raise WorkflowError(
            "'nextflow' must be in PATH in order to execute nextflow directive."
        )

    args = []
    if use_singularity:
        args.append("-profile singularity")
    elif use_conda:
        args.append("-profile conda")
    # TODO pass threads in case of single job
    # TODO limit parallelism in case of pipeline
    # TODO handle other resources

    add_parameter = lambda name, value: args.append("--{} {}".format(name, value))

    for name, files in input.items():
        if is_pipeline and any(f.is_remote for f in files):
            # TODO handle remote filesystems if input file uses remote provider
            raise WorkflowError(
                "Remote files are not yet supported in combination with nextflow pipelines."
            )
        # TODO how are multiple input files under a single arg usually passed to nextflow?
        add_parameter(name, files)
    for name, value in params.items():
        add_parameter(name, value)

    log = "2> {}".format(log) if log else ""

    cmd = "nextflow {args} {repo} {log}".format(args=" ".join(args), repo=repo, log=log)
    shell(cmd, bench_record=bench_record)
