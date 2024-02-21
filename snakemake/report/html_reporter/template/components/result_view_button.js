'use strict';

class ResultViewButton extends React.Component {
    render() {
        let result = results[this.props.resultPath];

        return this.getViewButton(this.props.resultPath, result);
    }

    getViewButton(resultPath, entry) {
        const mimeType = this.getResultMimeType(resultPath);
        let setView = this.props.app.setView;

        let props;

        switch (mimeType) {
            case "image/svg+xml":
            case "image/png":
            case "image/jpeg":
                props = {
                    href: "#",
                    onClick: function () {
                        setView({
                            content: "img",
                            contentPath: entry.data_uri
                        })
                    }
                };
                break;
            case "text/html":
                props = {
                    href: "#",
                    onClick: function () {
                        setView({
                            content: "html",
                            contentPath: entry.data_uri
                        })
                    }
                };
                break;
            case "application/pdf":
                props = {
                    href: "#",
                    onClick: function () {
                        setView({
                            content: "pdf",
                            contentPath: entry.data_uri
                        })
                    }
                };
                break;
            default:
                props = {
                    href: entry.data_uri,
                    download: entry.name,
                    target: "_blank"
                };
        }
        return e(
            Button,
            { iconName: "eye", ...props }
        );
    }

    getResultMimeType(resultPath) {
        return results[resultPath].mime_type
    }
}

ResultViewButton.propTypes = {
    resultPath: PropTypes.string.isRequired,
    app: PropTypes.object.isRequired
};