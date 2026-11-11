'use strict';

class VegaPlot extends React.Component {
    constructor(props) {
        super(props);
        this.container = React.createRef();
        this.embedded = null;
        this.spec = null;
        this.specWithoutLegends = null;
        this.hasLegend = false;
        this.baseWidth = null;
        this.baseHeight = null;
        this.scale = 100;
        this.showLegend = true;
        this.state = { error: null };
        this.showError = this.showError.bind(this);
    }

    componentDidMount() {
        const result = results[this.props.resultPath];
        const specPromise = result.vega_spec
            ? Promise.resolve(
                JSON.parse(LZString.decompressFromBase64(result.vega_spec))
            )
            : fetch(result.data_uri).then((response) => response.json());
        specPromise
            .then((spec) => {
                this.spec = spec;
                return this.renderPlot();
            })
            .catch(this.showError);
    }

    componentWillUnmount() {
        if (this.embedded && this.embedded.finalize) {
            this.embedded.finalize();
        }
    }

    showError(err) {
        this.setState({ error: "Could not render plot: " + err });
    }

    renderPlot() {
        const el = this.container.current;
        if (this.embedded && this.embedded.finalize) {
            this.embedded.finalize();
        }
        const options = { actions: true, forceActionsMenu: true };
        const embedPromise = this.showLegend
            ? vegaEmbed(el, this.spec, options)
            : vegaEmbed(el, this.specWithoutLegends, { ...options, mode: "vega" });
        return embedPromise.then((embedded) => {
            this.embedded = embedded;
            if (this.baseWidth === null) {
                // first rendering: remember the original size and look for legends
                this.baseWidth = embedded.view.width();
                this.baseHeight = embedded.view.height();
                const legends = embedded.vgSpec ? embedded.vgSpec.legends : null;
                this.hasLegend = Array.isArray(legends) && legends.length > 0;
                if (this.hasLegend) {
                    const copy = JSON.parse(JSON.stringify(embedded.vgSpec));
                    delete copy.legends;
                    this.specWithoutLegends = copy;
                }
            } else {
                this.applyScale();
            }
            this.addMenuControls(el);
        });
    }

    applyScale() {
        this.embedded.view
            .width(this.baseWidth * this.scale / 100)
            .height(this.baseHeight * this.scale / 100)
            .run();
    }

    addMenuControls(el) {
        const menu = el.querySelector(".vega-actions");
        if (!menu) {
            return;
        }
        menu.append(this.sizeRow());
        if (this.hasLegend) {
            menu.append(this.legendRow());
        }
    }

    menuRow() {
        const row = document.createElement("div");
        row.style.cssText =
            "padding: 8px 16px; font-family: sans-serif; font-size: 14px; " +
            "font-weight: 600; color: #434a56; white-space: nowrap; " +
            "display: flex; align-items: center; gap: 8px;";
        return row;
    }

    sizeRow() {
        const row = this.menuRow();
        const label = document.createElement("span");
        label.textContent = "Size";
        const slider = document.createElement("input");
        slider.type = "range";
        slider.min = 25;
        slider.max = 300;
        slider.step = 5;
        slider.value = this.scale;
        const value = document.createElement("span");
        value.textContent = this.scale + "%";

        // while dragging, only update the label
        slider.addEventListener("input", () => {
            value.textContent = slider.value + "%";
        });
        // re-render once, when the mouse button is released, and only if the size changed
        slider.addEventListener("change", () => {
            const scale = Number(slider.value);
            if (scale !== this.scale) {
                this.scale = scale;
                this.applyScale();
            }
        });

        row.append(label, slider, value);
        return row;
    }

    legendRow() {
        const row = this.menuRow();
        const label = document.createElement("label");
        label.style.cssText = "display: flex; align-items: center; gap: 8px; cursor: pointer;";
        const checkbox = document.createElement("input");
        checkbox.type = "checkbox";
        checkbox.checked = this.showLegend;
        checkbox.addEventListener("change", () => {
            this.showLegend = checkbox.checked;
            this.renderPlot().catch(this.showError);
        });
        label.append(checkbox, "Legend");
        row.append(label);
        return row;
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