def render_categories(results):
    return {
        cat.name: {
            subcat.name: [res.path for res in catresults]
            for subcat, catresults in subcats
        }
        for cat, subcats in results
    }
