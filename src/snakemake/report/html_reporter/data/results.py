import json
import os

from humanfriendly import format_size
from snakemake_interface_common.exceptions import WorkflowError

from snakemake.report.common import VEGA_MIME_TYPES
from snakemake.report.html_reporter.common import get_result_uri


def compress_vega_spec(res):
    """Return the lz-string compressed spec for Vega results, else None."""
    if res.mime not in VEGA_MIME_TYPES.values():
        return None
    try:
        import lzstring
    except ImportError:
        raise WorkflowError(
            "Python package lzstring must be installed to include Vega plots "
            "in reports."
        )
    with open(res.path, encoding="utf-8") as f:
        return lzstring.LZString().compressToBase64(f.read())


def get_data_uri(res, mode_embedded):
    # Embedded Vega specs are stored compressed (vega_spec), not as data URI.
    if mode_embedded and res.mime in VEGA_MIME_TYPES.values():
        return None
    return get_result_uri(res, mode_embedded)


def render_results(results, mode_embedded):
    return json.dumps(
        {
            str(res.path): {
                "name": res.name,
                "filename": res.filename,
                "labels": res.labels,
                "size": format_size(res.size),
                "caption": res.caption,
                "mime_type": res.mime,
                "job_properties": {
                    "rule": res.job.rule.name,
                    "wildcards": res.wildcards,
                    "params": res.params,
                },
                "vega_spec": compress_vega_spec(res),
                "data_uri": get_data_uri(res, mode_embedded),
            }
            for cat, subcats in results.items()
            for subcat, catresults in subcats.items()
            for res in catresults
        }
    )
