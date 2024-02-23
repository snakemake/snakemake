class AbstractResults extends React.Component {
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

    getLabels() {
        let first_index = {};
        this.getResults().map(function ([path, result]) {
            let i = 0;
            for (let key in result.labels) {
                if (!(key in first_index)) {
                    first_index[key] = i;
                }
                i += 1;
            }
        })
        let labels = Object.keys(first_index);

        return labels.sort(function (a, b) {
            return first_index[a] - first_index[b];
        });
    }

    isLabelled() {
        return this.getResults().every(function ([path, result]) {
            return result.labels;
        });
    }

    renderHeader() {
        if (this.isLabelled()) {
            return e(
                "tr",
                {},
                this.getLabels().map(function (label) {
                    return e(
                        "th",
                        { className: "text-left p-1 uppercase" },
                        label
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
        AbstractResults.propTypes = {
            app: PropTypes.object.isRequired,
        };

        let app = this.props.app;
        let labels;
        if (this.isLabelled()) {
            labels = this.getLabels();
        }
        return this.getResults().map(function ([path, entry]) {
            let actions = e(
                "td",
                { className: "p-1 text-right" },
                e(
                    "div",
                    { className: "inline-flex gap-1", role: "group" },
                    e(
                        ResultViewButton,
                        { resultPath: path, app: app }
                    ),
                    e(
                        Button,
                        {
                            href: "#",
                            onClick: () => app.showResultInfo(path),
                            iconName: "information-circle"
                        }
                    )
                )
            );

            let entryLabels;
            let key;
            if (labels !== undefined) {
                entryLabels = labels.map(function (label) {
                    return e(
                        "td",
                        { className: "p-1" },
                        entry.labels[label] || ""
                    );
                });
                key = labels.join();
            } else {
                entryLabels = e(
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
                    entryLabels,
                    actions
                )
            ];
        })
    }
}