'use strict';

class ListHeading extends React.Component {
    render() {
        return e(
            "li",
            { className: "uppercase font-bold p-1" },
            this.props.text
        );
    }
}