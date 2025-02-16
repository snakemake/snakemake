function arrayKey(array) {
    return array.join(",");
}

class AbstractResults extends React.Component {
    constructor(props) {
        super(props);
        this.state = { toggles: new Map(), data: this.getData() };
        this.toggleCallback = this.toggleCallback.bind(this);
    }

    render() {
        if (this.state.data.toggleLabels.size > 0) {
            this.initToggleStates(this.state.data.toggleLabels);
            return e(
                "div",
                {},
                e(
                    "div",
                    { className: "flex gap-2 rounded bg-slate-800" },
                    this.getToggleControls(this.state.data.toggleLabels),
                ),
                this.getResultsTable(this.state.data)
            )
        } else {
            return this.getResultsTable(this.state.data);
        }
    }

    initToggleStates(toggleLabels) {
        let toggles = this.state.toggles;
        toggles.clear();
        toggleLabels.forEach(function(value, key) {
            toggles.set(key, value[0]);
        })
    }

    getToggleControls(toggleLabels) {
        let toggleCallback = this.toggleCallback;
        return toggleLabels.entries().map(function(entry) {
            let [name, values] = entry;
            return e(
                Toggle,
                {
                    label: name,
                    values: values,
                    callback: function(selected) {
                        toggleCallback(name, selected);
                    }
                }
            )
        })
    }

    toggleCallback(name, selected) {
        let data = this.state.data;
        this.setState(function(prevState) {
            let toggles = new Map(prevState.toggles);
            toggles.set(name, selected);
            return { data: data, toggles };
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
        if (labels !== undefined) {
            labels.forEach(function (label) {
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

        let entries = new Map();
        let entryLabelValues = [];

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
            entryLabels: labels.filter((label) => toggleLabels.has(label)),
            entryLabelValues: entryLabelValues,
            toggleLabels: toggleLabels,
            entries: entries,
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
            let toggleLabels = Array.from(data.toggleLabels.keys().map((label) => state.toggles.get(label)));
            let entryPath = data.entries.get(arrayKey(entryLabels)).get(arrayKey(toggleLabels));
            console.log({
                toggleLabels,
                entryPath,
                entryLabels,
            });
            

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