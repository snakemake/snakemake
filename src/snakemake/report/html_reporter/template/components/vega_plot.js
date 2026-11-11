'use strict';

class VegaPlot extends React.Component {
    constructor(props) {
        super(props);
        this.container = React.createRef();
        this.embedded = null;
        this.spec = null;
        this.isVegaLite = false;
        this.hasLegend = false;
        this.baseWidth = null;
        this.baseHeight = null;
        this.scale = 100;
        this.showLegend = true;
        // "?vega_size=container" in the report URL makes container sizing the default
        this.fit = new URLSearchParams(window.location.search).get("vega_size") === "container";
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
                this.isVegaLite = /vega-lite/.test(spec.$schema || "");
                const legends = this.compiledSpec(spec).legends;
                this.hasLegend = Array.isArray(legends) && legends.length > 0;
                if (!this.isVegaLite) {
                    this.fit = false;
                }
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

    // plain Vega version of a spec (a copy)
    compiledSpec(spec) {
        return this.isVegaLite
            ? vegaLite.compile(spec).spec
            : JSON.parse(JSON.stringify(spec));
    }

    // spec and options for the current settings (legend shown or hidden, container size)
    embedArguments() {
        const spec = this.fit ? { ...this.spec, width: "container" } : this.spec;
        const options = {
            actions: { export: false, source: true, compiled: true, editor: true },
            forceActionsMenu: true
        };
        if (this.showLegend) {
            return [spec, options];
        }
        const vgSpec = this.compiledSpec(spec);
        delete vgSpec.legends;
        return [vgSpec, { ...options, mode: "vega" }];
    }

    renderPlot() {
        const el = this.container.current;
        if (this.embedded && this.embedded.finalize) {
            this.embedded.finalize();
        }
        // (the report centers its content in a flex box, so the wrapper must be widened too)
        el.style.display = this.fit ? "block" : "";
        el.style.width = this.fit ? "100%" : "";
        el.parentElement.style.width = this.fit ? "100%" : "";
        const [spec, options] = this.embedArguments();
        return vegaEmbed(el, spec, options).then((embedded) => {
            this.embedded = embedded;
            this.baseWidth = embedded.view.width();
            this.baseHeight = embedded.view.height();
            if (this.scale !== 100) {
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

    // the view to save: with container size and no manual resizing, the original size
    exportView() {
        if (this.fit && this.scale === 100) {
            const vgSpec = this.compiledSpec(this.spec);
            if (!this.showLegend) {
                delete vgSpec.legends;
            }
            const view = new vega.View(vega.parse(vgSpec), { renderer: "none" });
            return view.runAsync().then(() => view);
        }
        return Promise.resolve(this.embedded.view);
    }

    download(blob, extension) {
        const link = document.createElement("a");
        link.href = URL.createObjectURL(blob);
        link.download = "visualization." + extension;
        document.body.append(link);
        link.click();
        link.remove();
        setTimeout(() => URL.revokeObjectURL(link.href), 1000);
    }

    saveSvg() {
        return this.exportView()
            .then((view) => view.toSVG())
            .then((svg) => this.download(new Blob([svg], { type: "image/svg+xml" }), "svg"));
    }

    savePng() {
        return this.exportView()
            .then((view) => view.toCanvas(1))
            .then((canvas) => canvas.toBlob((blob) => this.download(blob, "png")));
    }

    addMenuControls(el) {
        const menu = el.querySelector(".vega-actions");
        if (!menu) {
            return;
        }
        const closeMenu = () => {
            const details = el.querySelector("details");
            if (details) {
                details.removeAttribute("open");
            }
        };
        const saveItem = (text, save) => {
            const item = document.createElement("a");
            item.textContent = text;
            item.addEventListener("click", () => {
                closeMenu();
                save().catch(this.showError);
            });
            return item;
        };
        menu.prepend(saveItem("Save as PNG", () => this.savePng()));
        menu.prepend(saveItem("Save as SVG", () => this.saveSvg()));
        menu.append(this.sizeRow());
        if (this.isVegaLite) {
            menu.append(this.fitRow());
        }
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

    checkboxRow(text, checked, onChange) {
        const row = this.menuRow();
        const label = document.createElement("label");
        label.style.cssText = "display: flex; align-items: center; gap: 8px; cursor: pointer;";
        const checkbox = document.createElement("input");
        checkbox.type = "checkbox";
        checkbox.checked = checked;
        checkbox.addEventListener("change", () => onChange(checkbox.checked));
        label.append(checkbox, text);
        row.append(label);
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

    fitRow() {
        return this.checkboxRow("Fit to container", this.fit, (checked) => {
            this.fit = checked;
            this.scale = 100;
            this.renderPlot().catch(this.showError);
        });
    }

    legendRow() {
        return this.checkboxRow("Legend", this.showLegend, (checked) => {
            this.showLegend = checked;
            this.renderPlot().catch(this.showError);
        });
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