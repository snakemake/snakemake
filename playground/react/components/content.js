'use strict';



class ContentDisplay extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        return e(
            "div",
            { className: "flex items-center justify-center min-h-screen z-0" },
            this.renderContent()
        )
    }

    renderContent() {
        switch (this.props.show) {
            case "rulegraph":
                return e(RuleGraph);
            case "stats":
                return e(Stats)
        }
    }
}