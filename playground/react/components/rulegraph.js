'use strict';

class RuleGraph extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        return e(
            "div",
            { id: "rulegraph" }
        )
    }

    componentDidMount() {
        let showRuleProperties = this.showRuleProperties;
        vegaEmbed("#rulegraph", rulegraph_spec).then(function (ret) {
            ret.view.addEventListener("click", function (event, item) {
                if (item && "rule" in item.datum) {
                    var rule = item.datum.rule;
                    showRuleProperties(rule);
                }
            });
        });
    }


    showRuleProperties() {

    }
}