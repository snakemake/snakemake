from abc import ABC, abstractmethod

class ForeignPipeline:
    def __init__(
        self, 
        input,
        params,
        parallelism,
        resources,
        log,
        use_singularity,
        use_conda,
        bench_record,
    ):
        self.input = input
        self.params = params
        self.parallelism = parallelism
        self.resources = resources
        self.log = log
        self.use_singularity = use_singularity
        self.use_conda = use_conda
        self.bench_record = bench_record

    @abstractmethod
    def run(self):
        return NotImplemented
