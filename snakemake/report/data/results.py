import json


def render_results(results):
    return json.dumps(
        {
            res.path: {
                "name": res.name,
                "filename": res.filename,
                "labels": res.labels,
                "size": res.size_mb,
                "caption": res.caption,
                "mime_type": res.mime,
                "job_properties": {
                    "rule": res.job.rule.name,
                    "wildcards": res.wildcards,
                    "params": res.params,
                },
                "data_uri": res.data_uri,
            }
            for cat, subcats in results.items()
            for subcat, catresults in subcats.items()
            for res in catresults
        }
    )
