'use strict';

class ListItem extends React.Component {
    render() {
        return e(
            "li",
            { className: "p-1", ...this.props },
            this.props.children
        );
    }
}