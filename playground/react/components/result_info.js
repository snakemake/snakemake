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

        if (result.columns !== undefined) {
            const columns = Object.keys(result.columns).sort();
            return [
                e(
                    ListHeading,
                    { text: "Result", key: "result" }
                ),
                e(
                    ListItem,
                    e(
                        "table",
                        {},
                        e(
                            "thead",
                            {},
                            e(
                                "tr",
                                {},
                                columns.map(function (column) {
                                    return e(
                                        "th",
                                        {},
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
                                        {},
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
                    { text: "Result", key: "result" }
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
        console.log(caption)
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
                        className: "p-1",
                        dangerouslySetInnerHTML: { __html: caption }
                    }
                )
            ];
        } else {
            return [];
        }
    }

    getRule() {
        return [
            e(
                ListHeading,
                { key: "rule-heading", text: "Snakemake rule" }
            ),
            e(
                ListItem,
                { key: "rulename" },
                this.getResult().job_properties.rule
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