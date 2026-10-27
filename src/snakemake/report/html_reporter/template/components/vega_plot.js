'use strict';

class VegaPlot extends React.Component {
    constructor(props) {
        super(props);
        this.container = React.createRef();
        this.view = null;
        this.baseWidth = null;
        this.baseHeight = null;
        this.state = { scale: 100, error: null };
        this.onScaleChange = this.onScaleChange.bind(this);
    }

    componentDidMount() {
        const el = this.container.current;
        const result = results[this.props.resultPath];
        const specPromise = result.vega_spec
            ? Promise.resolve(
                JSON.parse(LZString.decompressFromBase64(result.vega_spec))
            )
            : fetch(result.data_uri).then((response) => response.json());
        specPromise
            .then((spec) => vegaEmbed(el, spec, { actions: true }))
            .then((embedded) => {
                this.view = embedded.view;
                this.baseWidth = embedded.view.width();
                this.baseHeight = embedded.view.height();
            })
            .catch((err) => {
                this.setState({ error: "Could not render plot: " + err });
            });
    }

    onScaleChange(event) {
        const scale = Number(event.target.value);
        this.setState({ scale: scale });
        if (this.view) {
            this.view
                .width(this.baseWidth * scale / 100)
                .height(this.baseHeight * scale / 100)
                .run();
        }
    }

    render() {
        return e(
            "div",
            { className: "p-3" },
            e(
                "div",
                { className: "flex items-center gap-2 mb-3 text-sm" },
                e("span", {}, "Size"),
                e("input", {
                    type: "range",
                    min: 25,
                    max: 300,
                    step: 5,
                    value: this.state.scale,
                    onChange: this.onScaleChange
                }),
                e("span", {}, this.state.scale + "%")
            ),
            this.state.error ? e("div", {}, this.state.error) : null,
            e("div", { ref: this.container })
        );
    }
}

VegaPlot.propTypes = {
    resultPath: PropTypes.string.isRequired
};