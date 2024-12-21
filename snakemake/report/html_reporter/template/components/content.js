'use strict';



class ContentDisplay extends React.Component {
    render() {
        return e(
            "div",
            { className: "grow flex items-center justify-center h-screen overflow-auto" },
            this.renderContent()
        )
    }

    renderContent() {
        let setView = this.props.app.setView;
        ContentDisplay.propTypes = {
            app: PropTypes.shape({
                state: PropTypes.shape({
                    content: PropTypes.string.isRequired,
                    contentPath: PropTypes.string,
                    contentText: PropTypes.string
                }).isRequired,
                setView: PropTypes.func.isRequired
            }).isRequired
        };

        switch (this.props.app.state.content) {
            case "rulegraph":
                return e(
                    "div",
                    { className: "grow flex gap-3 p-3 items-start" },
                    e(
                        "div",
                        { className: "py-2" },
                        e(
                            "div",
                            {
                                className: "prose prose-sm max-w-lg",
                                dangerouslySetInnerHTML: { __html: workflow_desc }
                            }
                        ),
                        e(
                            "div",
                            { id: "brand" }
                        )
                    ),
                    e(RuleGraph, { setView: setView })
                );
            case "stats":
                return e(
                    "div",
                    { className: "p-3" },
                    e(Stats)
                );
            case "img":
                return e(
                    "div",
                    { className: "p-3" },
                    e(
                        "img",
                        { src: this.props.app.state.contentPath }
                    )
                );
            case "html":
            case "pdf":
                return e(
                    "iframe",
                    { src: this.props.app.state.contentPath, className: "w-full h-screen" }
                );
            case "text":
                return e(
                    "div",
                    { className: "p-3 w-full" },
                    e(
                        "pre",
                        { className: "whitespace-pre-line text-sm" },
                        this.props.app.state.contentText
                    )
                );
        }
    }
}