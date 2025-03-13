'use strict';

class ListItem extends React.Component {
    static propTypes = {
        children: PropTypes.node.isRequired,
    };

    render() {
        return e(
            "li",
            { className: "p-1", ...this.props },
            this.props.children
        );
    }
}