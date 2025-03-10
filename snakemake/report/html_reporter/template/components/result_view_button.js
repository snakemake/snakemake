'use strict';


class ButtonViewManager extends AbstractViewManager {
    constructor(app) {
        super();
        this.app = app;
    }

    handleImg(entry, resultPath) {
        let app = this.app;
        let props = {
            href: "#",
            onClick: function() {
                app.setView({
                    content: "img",
                    contentPath: entry.data_uri,
                    resultPath: resultPath
                });
            }
        };
        return this.renderButton(props);
    }

    handleHtml(entry, resultPath) {
        let app = this.app;
        let props = {
            href: "#",
            onClick: function() {
                app.setView({
                    content: "html",
                    contentPath: entry.data_uri,
                    resultPath: resultPath
                });
            }
        };
        return this.renderButton(props);
    }

    handlePdf(entry, resultPath) {
        let app = this.app;
        let props = {
            href: "#",
            onClick: function() {
                app.setView({
                    content: "pdf",
                    contentPath: entry.data_uri,
                    resultPath: resultPath
                });
            }
        };
        return this.renderButton(props);
    }

    handleDefault(entry, resultPath) {
        let props = {
            href: entry.data_uri,
            download: entry.name,
            target: "_blank"
        };
        return this.renderButton(props);
    }

    renderButton(props) {
        return e(
            Button,
            { iconName: "eye", ...props }
        );
    }
}


class ResultViewButton extends React.Component {
    constructor(props) {
        super(props);
        this.viewManager = new ButtonViewManager(this.props.app);
    }

    render() {
        return this.viewManager.handleSelectedResult(this.props.resultPath);
    }
}

ResultViewButton.propTypes = {
    resultPath: PropTypes.string.isRequired,
    app: PropTypes.object.isRequired
};