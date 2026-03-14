import json


def render_metadata(metadata):
    """Render the metadata to be displayed on the report landing page.

    Arguments
    ---------
    metadata: dict -- Metadata dictionary containing (user specified) metadata
    """
    return json.dumps(metadata)
