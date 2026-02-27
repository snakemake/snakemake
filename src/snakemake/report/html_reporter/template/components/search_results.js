'use strict';

class SearchResults extends AbstractResults {
    getSearchFunc() {
        let term = this.props.searchTerm;
        if (term.startsWith("re:")) {
            let regexp = new RegExp(term.slice(3));
            return function (text) {
                return regexp.test(text);
            }
        } else {
            return function (text) {
                return text.includes(term);
            }
        }
    }

    getResults() {
        let searchFunc = this.getSearchFunc();
        return Object.entries(results).map(function ([path, result]) {
            if (searchFunc(path)) {
                return [path, result];
            }
            const columns = result.columns || [];
            for (const columnValue in columns) {
                if (searchFunc(columnValue)) {
                    return [path, result];
                }
            }
            return undefined;
        }).filter(function (item) {
            return item !== undefined;
        })
    }

    getCategory() {
        return undefined;
    }

    getSubcategory() {
        return undefined;
    }

    getSearchTerm() {
        return this.props.searchTerm;
    }
}