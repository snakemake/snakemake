'use strict';

class Stats extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        return e(
            "div",
            { className: "flex" },
            e(
                "div",
                { id: "runtimes" }
            ),
            e(
                "div",
                { id: "timeline" }
            )
        )
    }

    componentDidMount() {
        vegaEmbed("#runtimes", runtimes_spec);
        vegaEmbed("#timeline", timeline_spec);
    }
}