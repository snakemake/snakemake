'use strict';

class AbstractResults extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        return e(
            "table",
            { className: "table-auto text-white text-sm w-full" },
            e(
                "thead",
                {},
                this.renderHeader(),
            ),
            e(
                "tbody",
                {},
                this.renderEntries()
            )
        )
    }

    getResults() {
        throw new Error("Unimplemented!");
    }

    getCategory() {
        throw new Error("Unimplemented!");
    }

    getSubcategory() {
        throw new Error("Unimplemented!");
    }

    getSearchTerm() {
        throw new Error("Unimplemented!");
    }

    getColumns() {
        return Array.from(new Set(this.getResults().map(function ([path, result]) { Object.keys(result.columns) }).flat())).sort();
    }

    isColumnBased() {
        return this.getResults().every(function ([path, result]) {
            "columns" in result
        });
    }

    renderHeader() {
        if (this.isColumnBased()) {
            return e(
                "tr",
                {},
                this.getColumns().map(function (column) {
                    return e(
                        "th",
                        { className: "text-left p-1" },
                        column
                    )
                }),
                e(
                    "th",
                    { className: "text-right p-1 w-fit" },
                )
            )
        } else {
            return e(
                "tr",
                {},
                e(
                    "th",
                    { className: "text-left p-1" },
                    "Filename"
                ),
                e(
                    "th",
                    { className: "text-right p-1 w-fit" },
                )
            )
        }
    }

    renderEntries() {
        let _this = this;
        let columns = undefined;
        if (this.isColumnBased()) {
            columns = this.getColumns();
        }
        return this.getResults().map(function ([path, entry]) {
            let actions = e(
                "td",
                { className: "p-1 text-right" },
                e(
                    "div",
                    { className: "inline-flex gap-1", role: "group" },
                    _this.getViewButton(path, entry),
                    _this.renderButton(
                        "information-circle",
                        {
                            href: "#",
                            onClick: () => _this.showResultInfo(path)
                        }
                    )
                )
            );

            let entryColumns = undefined;
            let key = undefined;
            if (columns !== undefined) {
                entryColumns = columns.map(function (column) {
                    return e(
                        "td",
                        { className: "p-1" },
                        entry.columns[column] || ""
                    );
                });
                key = columns.join();
            } else {
                entryColumns = e(
                    "td",
                    { className: "p-1" },
                    path
                );
                key = path;
            }

            return [
                e(
                    "tr",
                    { key: key },
                    entryColumns,
                    actions
                )
            ];
        })
    }

    getViewButton(resultPath, entry) {
        const mimeType = this.getResultMimeType(resultPath);
        let setView = this.props.setView;

        let props = undefined;

        switch (mimeType) {
            case "image/svg+xml":
            case "image/png":
            case "image/jpeg":
                props = {
                    href: "#",
                    onClick: function () {
                        setView({
                            content: "img",
                            contentPath: entry.data_uri
                        })
                    }
                };
                break;
            case "text/html":
                props = {
                    href: "#",
                    onClick: function () {
                        setView({
                            content: "html",
                            contentPath: entry.data_uri
                        })
                    }
                };
                break;
            default:
                props = {
                    href: entry.data_uri,
                    download: entry.name
                };
        }
        return this.renderButton("eye", props);
    }

    renderButton(iconName, props) {
        return e(
            "a",
            { type: "button", className: `transition-all inline-block p-1 text-emerald-500 rounded hover:bg-slate-800`, ...props },
            e(Icon, { iconName: iconName })
        )
    }

    showResultInfo(resultPath) {
        this.props.setView({ navbarMode: "resultinfo", resultPath: resultPath, category: this.getCategory(), subcategory: this.getSubcategory(), searchTerm: this.getSearchTerm() });
    }

    getResultMimeType(resultPath) {
        return results[resultPath].mime_type
    }
}