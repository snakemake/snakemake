'use strict';

class ResultInfo extends React.Component {
    render() {
        return e(
            "ul",
            {},
            this.getDescriptor(),
            this.getCaption(),
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
                    { text: "Result" }
                ),
                e(
                    "li",
                    {},
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
                    { text: "Result" }
                ),
                e(
                    "li",
                    { className: "p-1" },
                    this.props.resultPath
                )
            ];
        }
    }

    getCaption() {
        const caption = this.getResult().caption;
        if (caption !== undefined) {
            return [
                e(
                    ListHeading,
                    { text: "Description" }
                ),
                e(
                    "li",
                    {
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
                { text: "Snakemake rule" }
            ),
            e(
                ListItem,
                {},
                this.getResult().job_properties.rule
            )
        ]
    }

    getParamsAndWildcards() {
        const jobProperties = this.getResult().job_properties;
        console.log(jobProperties)
        if (jobProperties.params || jobProperties.wildcards) {
            let items = [
                e(
                    ListHeading,
                    { text: "Job parameters" }
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
}