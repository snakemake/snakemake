import json

import snakemake


def render_metadata(metadata):
    """Render the metadata from to be displayed on the report landing page.

    Arguments
    ---------
    metadata: MetadataRecord -- Metadata object containing (user specified) metadata
    """
    result_dict = {"Snakemake version": metadata.snakemake_version}
    if metadata.metadata_dict:
        result_dict.update(metadata.metadata_dict)
    return json.dumps(result_dict)