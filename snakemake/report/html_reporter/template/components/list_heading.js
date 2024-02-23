'use strict';

class ListHeading extends React.Component {
    static propTypes = {
        text: PropTypes.string.isRequired
    };

    render() {
        return e(
            "li",
            { className: "uppercase font-bold p-1" },
            this.props.text
        );
    }
}