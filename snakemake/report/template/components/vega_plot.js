'use strict';

class Stats extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        return e(
            "div",
            { id: "vega-plot" }
        )
    }

    componentDidMount() {
        vegaEmbed("#vega-plot", this.props.app.state.contentPath);
    }
}