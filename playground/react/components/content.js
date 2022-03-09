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
        let setView = this.props.app.setView;
        switch (this.props.app.state.content) {
            case "rulegraph":
                return e(RuleGraph, { setView: setView });
            case "stats":
                return e(Stats)
        }
    }
}