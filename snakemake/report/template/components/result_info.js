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

        if (result.columns) {
            const columns = Object.keys(result.columns).sort();
            return [
                e(
                    ListItem,
                    {},
                    e(
                        "table",
                        { className: "table-auto text-white text-sm" },
                        e(
                            "thead",
                            {},
                            e(
                                "tr",
                                {},
                                columns.map(function (column) {
                                    return e(
                                        "th",
                                        { className: "text-left uppercase pr-2" },
                                        column
                                    );
                                })
                            )
                        ),
                        e(
                            "tbody",
                            {},
                            e(
                                "tr",
                                {},
                                columns.map(function (column) {
                                    const value = result.columns[column];
                                    return e(
                                        "td",
                                        { className: "pr-2" },
                                        value
                                    );
                                })
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
                    this.props.resultPath
                )
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
        const setView = this.props.setView;
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