import json
import os

from humanfriendly import format_size

from snakemake.report.html_reporter.common import get_result_uri


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
                "data_uri": get_result_uri(res, mode_embedded),
            }
            for cat, subcats in results.items()
            for subcat, catresults in subcats.items()
            for res in catresults
        }
    )
