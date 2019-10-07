__author__ = "Alban Gaignard, Johannes Köster"
__copyright__ = "Copyright 2019, Alban Gaignard, Johannes Köster"
__license__ = "MIT"

import itertools

from snakemake.io import is_flagged

def create_research_object(dag, path):
    """Create a research object from workflow persistence information."""

    persistence = dag.workflow.persistence
    recorded_files = set()
    for job in dag.jobs:
        for f in itertools.chain(job.expanded_output, job.input):
            if f in recorded_files:
                continue
            
            # include metadata into the research object
            meta = persistence.metadata(f)

            if is_flagged(f, "report"):
                # include content into the research object
                