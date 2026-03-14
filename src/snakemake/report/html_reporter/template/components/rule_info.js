'use strict';

function flattenList(inputList) {
    const flattened = [];

    for (const item of inputList) {
        if (typeof item === 'object' && !Array.isArray(item)) {
            for (const [key, value] of Object.entries(item)) {
                if (Array.isArray(value)) {
                    flattened.push(...value.map(subItem => `${key}: ${subItem}`));
                } else {
                    flattened.push(`${key}: ${value}`);
                }
            }
        } else {
            flattened.push(item);
        }
    }

    return flattened;
}

class RuleInfo extends React.Component {
    static propTypes = {
        rule: PropTypes.object.isRequired,
    };

    render() {
        let rule = rules[this.props.rule];
        if (rule === undefined) {
            return e(
                "span",
                { className: "p-1" },
                `No metadata available for rule ${this.props.rule}`
            );
        }

        return e(
            "ol",
            {},
            this.renderItems("Input", rule.input, {}, false),
            this.renderItems("Output", rule.output),
            this.renderSoftware(),
            this.renderItems("Container", [rule.container_img_url]),
            this.renderCode(),
        )
    }

    renderSoftware() {
        let rule = rules[this.props.rule];
        if (rule.conda_env) {
            return this.renderItems("Software", flattenList(rule.conda_env.dependencies));
        } else {
            return [];
        }
    }

    renderCode() {
        let rule = rules[this.props.rule];
        if (rule.code.length) {
            return [
                e(
                    ListHeading,
                    { key: "code-heading", text: "Code" }
                ),
                e(
                    ListItem,
                    {
                        key: "code",
                        className: "p-1",
                        dangerouslySetInnerHTML: { __html: rule.code }
                    }
                )
            ];
        } else {
            return [];
        }
    }

    renderItems(heading, items, props = {}, margin = true) {
        if (items.length && items.every((item) => item !== undefined)) {
            let headingProps = {};
            if (margin) {
                headingProps = { className: "" }
            }
            return [
                e(
                    ListHeading,
                    { text: heading, ...headingProps }
                ),
                items.map(function (item) {
                    return e(
                        ListItem,
                        props,
                        item
                    );
                })
            ];
        } else {
            return [];
        }
    }
}