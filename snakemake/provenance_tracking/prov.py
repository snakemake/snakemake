from snakemake.provenance_tracking import ProvenanceTracker

class PROV(ProvenanceTracker):
    def record_metadata(self, output_path, metadata):
        """store metadata for given output in PROV"""
        pass

    def record_content(self, output_path):
        # do not store the file content in PROV
        pass