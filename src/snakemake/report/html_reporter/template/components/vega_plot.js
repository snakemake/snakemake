'use strict';

class VegaPlot extends React.Component {
    constructor(props) {
        super(props);
        this.container = React.createRef();
    }

    componentDidMount() {
        const el = this.container.current;
        fetch(this.props.specUri)
            .then((response) => response.json())
            .then((spec) => vegaEmbed(el, spec, { actions: true }))
            .catch((err) => {
                el.textContent = "Could not render plot: " + err;
            });
    }

    render() {
        return e("div", { className: "p-3", ref: this.container });
    }
}

VegaPlot.propTypes = {
    specUri: PropTypes.string.isRequired
};