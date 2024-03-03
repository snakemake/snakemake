import base64
import mimetypes
import os
from pathlib import Path
from snakemake.common import is_local_file
from snakemake.logging import logger

from snakemake_interface_common.exceptions import WorkflowError


def data_uri(data, filename, encoding="utf8", mime="text/plain"):
    """Craft a base64 data URI from file with proper encoding and mimetype."""
    data = base64.b64encode(data)
    return f'data:{mime};charset={encoding};filename={filename};base64,{data.decode("utf-8")}'


def mime_from_file(file):
    mime, encoding = mimetypes.guess_type(file)
    if mime is None:
        mime = "text/plain"
        logger.info(
            "Could not detect mimetype for {}, assuming text/plain.".format(file)
        )
    return mime, encoding


def data_uri_from_file(file, defaultenc="utf8"):
    """Craft a base64 data URI from file with proper encoding and mimetype."""
    if isinstance(file, Path):
        file = str(file)
    mime, encoding = mime_from_file(file)
    if encoding is None:
        encoding = defaultenc
    with open(file, "rb") as f:
        return data_uri(f.read(), os.path.basename(file), encoding, mime)
