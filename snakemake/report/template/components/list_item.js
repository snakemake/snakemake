'use strict';

class ListItem extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        return e(
            "li",
            { className: "p-1", ...this.props },
            this.props.children
        );
    }
}