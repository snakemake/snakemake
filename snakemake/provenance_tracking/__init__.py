__author__ = "Alban Gaignard, Johannes Köster"
__copyright__ = "Copyright 2019, Alban Gaignard, Johannes Köster"
__license__ = "MIT"

import itertools
from abc import ABC, abstractmethod

from snakemake.io import is_flagged


class ProvenanceTracker:
    @abstractmethod
    def record_metadata(self, output_path, metadata):
        pass

    @abstractmethod
    def record_content(self, output_path):
        pass

    def collect(self):
        """Collect information from workflow persistence storage."""

        persistence = dag.workflow.persistence
        recorded_files = set()
        for job in dag.jobs:
            for f in itertools.chain(job.expanded_output, job.input):
                if f in recorded_files:
                    print(f"RECORDED file {f}")
                    continue

                meta = persistence.metadata(f)

                # include metadata
                self.record_metadata(f, meta)

                if is_flagged(f, "report"):
                    # include content
                    self.record_content(f)
