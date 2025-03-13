const e = React.createElement;

function isNoResults() {
    return Object.keys(categories).length == 0;
}

function isSingleCategory() {
    return Object.keys(categories).length == 1;
}

function isSingleDefaultCategory() {
    return isSingleCategory() && Object.keys(categories)[0] == "Other";
}

function isSingleSubcategory(category) {
    return Object.keys(categories[category]).length == 1;
}

function getResultMimeType(resultPath) {
    return results[resultPath].mime_type
}