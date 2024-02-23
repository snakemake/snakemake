import json


def render_categories(results):
    return json.dumps(
        {
            cat.name: {
                subcat.name: [str(res.path) for res in catresults]
                for subcat, catresults in subcats.items()
                if catresults
            }
            for cat, subcats in results.items()
            if subcats and sum(len(catresults) for catresults in subcats.values())
        }
    )
