'use strict';

class ResultInfo extends React.Component {
    render() {
        return e(
            "ul",
            {},
            this.getDescriptor(),
            this.getCaption(),
            this.getSize(),
            this.getRule(),
            this.getParamsAndWildcards(),
        )
    }

    getResult() {
        return results[this.props.resultPath];
    }

    getDescriptor() {
        let result = this.getResult();
        let resultPath = this.props.resultPath;
        let app = this.props.app;

        if (result.labels) {
            const labels = Object.keys(result.labels).sort((a, b) => a.localeCompare(b));
            return [
                e(
                    ListItem,
                    {},
                    e(
                        "table",
                        { className: "table-auto text-white text-sm items-center" },
                        e(
                            "thead",
                            {},
                            e(
                                "tr",
                                {},
                                labels.map(function (label) {
                                    return e(
                                        "th",
                                        { className: "text-left uppercase pr-2 whitespace-nowrap" },
                                        label
                                    );
                                }),
                                e("th", {})
                            )
                        ),
                        e(
                            "tbody",
                            {},
                            e(
                                "tr",
                                {},
                                labels.map(function (label, index) {
                                    const value = result.labels[label];
                                    let item = value;
                                    if (index == labels.length - 1) {
                                        item = e(
                                            "span",
                                            { className: "flex items-center gap-2" },
                                            e("span", {}, value),
                                            e(
                                                ResultViewButton,
                                                { resultPath: resultPath, app: app }
                                            )
                                        );
                                    }

                                    return e(
                                        "td",
                                        { className: "pr-2" },
                                        item
                                    );
                                }),
                            )
                        )
                    )
                )
            ];
        } else {
            return [
                e(
                    ListHeading,
                    { text: "Path", key: "result" }
                ),
                e(
                    ListItem,
                    { key: "path" },
                    e(
                        "span",
                        { className: "flex items-center gap-2" },
                        e(
                            "span",
                            {},
                            this.props.resultPath
                        ),
                        e(
                            ResultViewButton,
                            { resultPath: this.props.resultPath, app: app }
                        )
                    )
                ),
            ];
        }
    }

    getCaption() {
        const caption = this.getResult().caption;
        if (caption) {
            return [
                e(
                    ListHeading,
                    { key: "caption-heading", text: "Description" }
                ),
                e(
                    ListItem,
                    {
                        key: "caption",
                        className: "p-1 prose prose-invert prose-sm",
                        dangerouslySetInnerHTML: { __html: caption }
                    }
                )
            ];
        } else {
            return [];
        }
    }

    getRule() {
        const setView = this.props.app.setView;
        const rule = this.getResult().job_properties.rule;
        return [
            e(
                ListHeading,
                { key: "rule-heading", text: "Snakemake rule" }
            ),
            e(
                ListItem,
                { key: "rulename" },
                e(
                    "span",
                    { className: "flex items-center gap-1" },
                    [
                        e(
                            "span",
                            {},
                            rule
                        ),
                        e(
                            "a",
                            {
                                type: "button",
                                href: "#",
                                className: `transition-all inline-block p-1 text-emerald-500 rounded hover:bg-slate-800`,
                                onClick: () => setView({ navbarMode: "ruleinfo", ruleinfo: rule })
                            },
                            e(Icon, { iconName: "information-circle" })
                        )
                    ]
                )
            )
        ]
    }

    getParamsAndWildcards() {
        const jobProperties = this.getResult().job_properties;
        if (jobProperties.params || jobProperties.wildcards) {
            let items = [
                e(
                    ListHeading,
                    { key: "params-heading", text: "Job parameters" }
                ),
            ];
            if (jobProperties.wildcards) {
                items.push(e(
                    ListItem,
                    { key: "wildcards" },
                    jobProperties.wildcards
                ));
            }
            if (jobProperties.params) {
                items.push(e(
                    ListItem,
                    { key: "params" },
                    jobProperties.params
                ));
            }
            return items;
        } else {
            return [];
        }
    }

    getSize() {
        return [
            e(
                ListHeading,
                { key: "size-heading", text: "Size" }
            ),
            e(
                ListItem,
                { key: "size" },
                this.getResult().size
            )
        ];
    }
}

ResultInfo.propTypes = {
    resultPath: PropTypes.string.isRequired,
    app: PropTypes.shape({
        state: PropTypes.shape({
            ruleinfo: PropTypes.string,
            subcategory: PropTypes.string,
            category: PropTypes.string,
            searchTerm: PropTypes.string,
            resultPath: PropTypes.string,
            hideNavbar: PropTypes.bool.isRequired,
            navbarMode: PropTypes.string.isRequired
        }).isRequired,
        setView: PropTypes.func.isRequired
    }).isRequired
};