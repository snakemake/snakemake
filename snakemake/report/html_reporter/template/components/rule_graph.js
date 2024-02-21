'use strict';

class RuleGraph extends React.Component {
    constructor(props) {
        super(props);
        this.showRuleProperties = this.showRuleProperties.bind(this);
    }

    render() {
        return e(
            "div",
            { id: "rulegraph", className: "max-h-screen py-2" }
        )
    }

    componentDidMount() {
        let showRuleProperties = this.showRuleProperties;
        vegaEmbed("#rulegraph", rulegraph_spec).then(function (ret) {
            ret.view.addEventListener("click", function (event, item) {
                if (item && "rule" in item.datum) {
                    let rule = item.datum.rule;
                    showRuleProperties(rule);
                }
            });
        });
    }

    showRuleProperties(rule) {
        this.props.setView({
            navbarMode: "ruleinfo",
            ruleinfo: rule
        });
    }
}

RuleGraph.propTypes = {
    setView: PropTypes.func.isRequired
};