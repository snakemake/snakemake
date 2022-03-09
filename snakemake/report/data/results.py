def render_results(results):
    return {
        res.path: {
            "name": res.name,
            "filename": res.filename,
            "size": res.size,
            "caption": res.caption,
            "job_properties": {
                "rule": res.rule,
                "wildcards": res.wildcards,
                "params": res.params,
            },
            "data_uri": res.data_uri,
        }
        for cat, subcats in results.items()
        for subcat, catresults in subcats
        for res in catresults
    }
