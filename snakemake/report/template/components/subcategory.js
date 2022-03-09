'use strict';

class Subcategory extends AbstractResults {
    getResults() {
        return categories[this.props.category][this.props.subcategory].map(function (path) {
            return [path, results[path]];
        });
    }

    getCategory() {
        return this.props.category;
    }

    getSubcategory() {
        return this.props.subcategory;
    }

    getSearchTerm() {
        return undefined;
    }
}