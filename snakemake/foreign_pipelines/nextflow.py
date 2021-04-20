import shutil
from itertools import chain
from snakemake.foreign_pipelines import ForeignPipeline

from snakemake.exceptions import WorkflowError
from snakemake.shell import shell
from snakemake.foreign_pipelines import ForeignPipeline


class NextflowPipeline(ForeignPipeline):
    def __init__(
        self,
        repo,
        input,
        params,
        parallelism,
        resources,
        log,
        use_singularity,
        use_conda,
        bench_record,
    ):
        super().__init__(input, params, parallelism, resources, log, use_singularity, use_conda, bench_record)
        self.repo = repo
    
    def run(self):
        if shutil.which("nextflow") is None:
            raise WorkflowError(
                "'nextflow' must be in PATH in order to execute nextflow directive."
            )
        
        args = []
        if self.use_singularity:
            args.append("-profile singularity")
        elif self.use_conda:
            args.append("-profile conda")
        # TODO pass threads in case of single job
        # TODO limit parallelism in case of pipeline
        # TODO handle other resources

        add_parameter = lambda name, value: args.append("--{} {}".format(name, value))

        for name, files in input.items():
            if any(f.is_remote for f in files):
                # TODO handle remote filesystems if input file uses remote provider
                raise WorkflowError(
                    "Remote files are not yet supported in combination with nextflow pipelines."
                )
            # TODO how are multiple input files under a single arg usually passed to nextflow?
            add_parameter(name, files)
        for name, value in self.params.items():
            add_parameter(name, value)

        log = "2> {}".format(log) if log else ""

        cmd = "nextflow {args} {repo} {log}".format(args=" ".join(args), repo=self.repo, log=log)
        shell(cmd, bench_record=self.bench_record)




