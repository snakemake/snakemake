function arrayKey(array) {
    return array.join(",");
}

class ToggleViewManager extends AbstractViewManager {
    constructor(app) {
        super();
        this.app = app;
    }

    handleImg(entry, resultPath) {
        this.app.setView({
            content: "img",
            contentPath: entry.data_uri,
            resultPath: resultPath
        });
    }

    handleHtml(entry, resultPath) {
        this.app.setView({
            content: "html",
            contentPath: entry.data_uri,
            resultPath: resultPath
        });
    }

    handlePdf(entry, resultPath) {
        this.app.setView({
            content: "pdf",
            contentPath: entry.data_uri,
            resultPath: resultPath
        });
    }

    handleDefault(entry, resultPath) {
        // do nothing in this case
    }
}

class AbstractResults extends React.Component {
    constructor(props) {
        super(props);
        let data = this.getData();
        let toggles = this.getInitToggleState(data.toggleLabels);
        this.state = { toggles, data };
        this.toggleCallback = this.toggleCallback.bind(this);
        this.toggleViewManager = new ToggleViewManager(this.props.app);
    }

    render() {
        if (this.state.data.toggleLabels.size > 0) {
            return e(
                "div",
                {},
                e(
                    "div",
                    { className: "p-2 flex flex-wrap gap-2 rounded bg-slate-800 text-xs" },
                    this.getToggleControls(this.state.data.toggleLabels),
                ),
                this.getResultsTable(this.state.data)
            )
        } else {
            return this.getResultsTable(this.state.data);
        }
    }

    getInitToggleState(toggleLabels) {
        let toggles = new Map();
        toggleLabels.forEach(function (value, key) {
            // Prefer "yes" as initial value if present, otherwise use the first value.
            let initialValue = value[1] === "yes" ? value[1] : value[0];
            toggles.set(key, initialValue);
        })
        return toggles;
    }

    getToggleControls(toggleLabels) {
        let toggleCallback = this.toggleCallback;
        let toggleState = this.state.toggles;
        return toggleLabels.entries().map(function (entry) {
            let [name, values] = entry;
            return e(
                Toggle,
                {
                    label: name,
                    values: values,
                    defaultValue: toggleState.get(name),
                    callback: function (selected) {
                        toggleCallback(name, selected);
                    }
                }
            )
        })
    }

    toggleCallback(name, selected) {
        let data = this.state.data;
        let _this = this;
        this.setState(function (prevState) {
            let toggles = new Map(prevState.toggles);
            toggles.set(name, selected);
            return { data: data, toggles };
        }, function () {
            if (_this.state.data.resultPathsToEntryLabels.has(_this.props.app.state.resultPath)) {
                let toggleLabels = Array.from(data.toggleLabels.keys()).map((label) => _this.state.toggles.get(label));
                let entryLabels = _this.state.data.resultPathsToEntryLabels.get(_this.props.app.state.resultPath);
                let targetPath = _this.state.data.entries.get(arrayKey(entryLabels)).get(arrayKey(toggleLabels));
                _this.toggleViewManager.handleSelectedResult(targetPath);
            }
        });
    }

    getResultsTable(data) {
        return e(
            "table",
            { className: "table-auto text-white text-sm w-full" },
            e(
                "thead",
                {},
                this.renderHeader(data),
            ),
            e(
                "tbody",
                {},
                this.renderEntries(data),
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

    getData() {
        let labels;
        if (this.isLabelled()) {
            labels = this.getLabels();
        }

        let results = this.getResults();

        let toggleLabels = new Map();
        // If there are at least two labels, consider all but the first label for
        // being shown as toggles. The latter is possible if a label
        // has exactly two values, each of which occur in half of the results.
        // Example: a plot which is created twice for each sample, once with and 
        // once without legend (for inclusion in larger panel figures where a 
        // repeating legend would be superfluous).
        if (labels !== undefined && labels.length > 1) {
            labels.slice(1).forEach(function (label) {
                let values = results.map(function ([path, entry]) {
                    return entry.labels[label];
                }).filter(function (value) {
                    return value !== undefined;
                });

                let uniqueValues = [...new Set(values)];
                if (
                    uniqueValues.length == 2 &&
                    uniqueValues.every(value => values.filter(v => v === value).length === results.length / 2)
                ) {
                    uniqueValues.sort();
                    toggleLabels.set(label, uniqueValues);
                }
            });
        }
        // Only allow one toggle label for now, in order to avoid confusion in the UI
        if (toggleLabels.size > 1) {
            toggleLabels = new Map();
        }

        let entries = new Map();
        let entryLabelValues = [];
        let resultPathsToEntryLabels = new Map();

        results.forEach(function ([path, entry]) {
            let entryLabels = [];
            let entryToggleLabels = [];
            if (labels !== undefined) {
                labels.forEach(function (label) {
                    if (!toggleLabels.has(label)) {
                        entryLabels.push(entry.labels[label] || "");
                    } else {
                        entryToggleLabels.push(entry.labels[label]);
                    }
                });
            } else {
                entryLabels = [path];
            }

            let key = arrayKey(entryLabels);

            if (!entries.has(key)) {
                entries.set(key, new Map());
                entryLabelValues.push(entryLabels);
            }

            entries.get(key).set(arrayKey(entryToggleLabels), path);
            resultPathsToEntryLabels.set(path, entryLabels);
        });

        entryLabelValues = entryLabelValues.sort(function (aLabels, bLabels) {
            // sort labels lexicographically, first element is the most important
            for (let i = 0; i < aLabels.length; i++) {
                let comparison = aLabels[i].localeCompare(bLabels[i]);
                if (comparison !== 0) {
                    return comparison;
                }
            }
            return 0;
        });

        if (labels === undefined) {
            labels = ["File"];
        }

        return {
            entryLabels: labels.filter((label) => !toggleLabels.has(label)),
            entryLabelValues,
            toggleLabels,
            entries,
            resultPathsToEntryLabels
        }
    }

    isLabelled() {
        return this.getResults().every(function ([path, result]) {
            return result.labels;
        });
    }

    renderHeader(data) {
        return e(
            "tr",
            {},
            data.entryLabels.map(function (label) {
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
    }

    renderEntries(data) {
        AbstractResults.propTypes = {
            app: PropTypes.object.isRequired,
        };

        let app = this.props.app;
        let state = this.state;

        return data.entryLabelValues.map(function (entryLabels) {
            let toggleLabels = Array.from(data.toggleLabels.keys()).map((label) => state.toggles.get(label));
            let entryPath = data.entries.get(arrayKey(entryLabels)).get(arrayKey(toggleLabels));

            let actions = e(
                "td",
                { className: "p-1 text-right" },
                e(
                    "div",
                    { className: "inline-flex gap-1", role: "group" },
                    e(
                        ResultViewButton,
                        { resultPath: entryPath, app: app }
                    ),
                    e(
                        Button,
                        {
                            href: "#",
                            onClick: () => app.showResultInfo(entryPath),
                            iconName: "information-circle"
                        }
                    )
                )
            );

            return [
                e(
                    "tr",
                    { key: entryLabels.join(",") },
                    entryLabels.map(function (labelValue) {
                        return e(
                            "td",
                            { className: "p-1" },
                            labelValue
                        );
                    }),
                    actions
                )
            ];
        });
    }
}
