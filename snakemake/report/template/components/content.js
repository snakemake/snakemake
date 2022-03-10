'use strict';



class ContentDisplay extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        return e(
            "div",
            { className: "flex items-center justify-center min-h-screen z-0 p-3" },
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
            case "img":
                return e(
                    "img",
                    { src: this.props.app.state.contentPath }
                )
            case "html":
                return e(
                    "iframe",
                    { src: this.props.add.state.contentPath, className: "w-screen h-screen" }
                )
        }
    }
}