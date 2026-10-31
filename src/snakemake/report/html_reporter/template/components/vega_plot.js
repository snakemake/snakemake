'use strict';

class VegaPlot extends React.Component {
    constructor(props) {
        super(props);
        this.container = React.createRef();
        this.state = { error: null };
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
            .then((spec) => vegaEmbed(el, spec, { actions: true, forceActionsMenu: true }))
            .then((embedded) => this.addSizeSlider(el, embedded.view))
            .catch((err) => {
                this.setState({ error: "Could not render plot: " + err });
            });
    }

    addSizeSlider(el, view) {
        const menu = el.querySelector(".vega-actions");
        if (!menu) {
            return;
        }
        const baseWidth = view.width();
        const baseHeight = view.height();

        const row = document.createElement("div");
        row.style.cssText =
            "padding: 8px 16px; font-family: sans-serif; font-size: 14px; " +
            "font-weight: 600; color: #434a56; white-space: nowrap; " +
            "display: flex; align-items: center; gap: 8px;";
        const label = document.createElement("span");
        label.textContent = "Size";
        const slider = document.createElement("input");
        slider.type = "range";
        slider.min = 25;
        slider.max = 300;
        slider.step = 5;
        slider.value = 100;
        const value = document.createElement("span");
        value.textContent = "100%";

        slider.addEventListener("input", function () {
            const scale = Number(slider.value);
            value.textContent = scale + "%";
            view.width(baseWidth * scale / 100)
                .height(baseHeight * scale / 100)
                .run();
        });
        slider.addEventListener("pointerdown", function () {
            menu.style.opacity = "0.3";
        });
        document.addEventListener("pointerup", function () {
            menu.style.opacity = "";
        });
        row.append(label, slider, value);
        menu.append(row);
    }

    render() {
        return e(
            "div",
            { className: "p-3" },
            this.state.error ? e("div", {}, this.state.error) : null,
            e("div", { ref: this.container })
        );
    }
}

VegaPlot.propTypes = {
    resultPath: PropTypes.string.isRequired
};