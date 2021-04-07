__author__ = "Alban Gaignard, Johannes Köster"
__copyright__ = "Copyright 2019, Alban Gaignard, Johannes Köster"
__license__ = "MIT"

from snakemake.logging import logger
from snakemake.provenance_tracking import ProvenanceTracker


class PROV(ProvenanceTracker):
    def record_metadata(self, output_path, metadata):
        """store metadata for given output in PROV"""
        logger.info(f"RECORD METADATA {metadata} to {output_path}")

    def record_content(self, output_path):
        # do not store the file content in PROV
        logger.info(f"RECORD METADATA to {output_path}")
