import json


def render_metadata(metadata):
    """Render the metadata from to be displayed on the report landing page.

    Arguments
    ---------
    metadata: MetadataRecord -- Metadata object containing (user specified) metadata
    """
    if metadata.metadata_dict:
        return json.dumps(metadata.metadata_dict)
    else:
        return json.dumps({})
