'use strict';

class RuleInfo extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        let rule = rules[this.props.rule];
        return e(
            "ol",
            {},
            this.renderItems("Input", rule.input, {}, false),
            this.renderItems("Output", rule.input),
            this.renderItems("Software", rule.conda),
            this.renderItems("Container", [rule.container]),
            this.renderItems("Code", rule.code),
        )
    }

    renderItems(heading, items, props = {}, margin = true) {
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
    }
}