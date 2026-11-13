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
        // by default, the plot is scaled up to use the available space
        this.fit = true;
        this.fitCheckbox = null;
        this.resizeFrame = null;
        this.state = { error: null };
        this.showError = this.showError.bind(this);
        this.onWindowResize = this.onWindowResize.bind(this);
    }

    componentDidMount() {
        window.addEventListener("resize", this.onWindowResize);
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
                return this.renderPlot();
            })
            .catch(this.showError);
    }

    componentWillUnmount() {
        window.removeEventListener("resize", this.onWindowResize);
        if (this.resizeFrame !== null) {
            cancelAnimationFrame(this.resizeFrame);
        }
        if (this.embedded && this.embedded.finalize) {
            this.embedded.finalize();
        }
    }

    showError(err) {
        this.setState({ error: "Could not render plot: " + err });
    }

    onWindowResize() {
        if (!this.fit || !this.embedded || this.resizeFrame !== null) {
            return;
        }
        this.resizeFrame = requestAnimationFrame(() => {
            this.resizeFrame = null;
            this.fitToWindow().catch(this.showError);
        });
    }

    // plain Vega version of a spec (a copy)
    compiledSpec(spec) {
        return this.isVegaLite
            ? vegaLite.compile(spec).spec
            : JSON.parse(JSON.stringify(spec));
    }

    // spec and options for the current settings (legend shown or hidden)
    embedArguments() {
        const options = {
            actions: { export: false, source: true, compiled: true, editor: true },
            forceActionsMenu: true
        };
        if (this.showLegend) {
            return [this.spec, options];
        }
        const vgSpec = this.compiledSpec(this.spec);
        delete vgSpec.legends;
        return [vgSpec, { ...options, mode: "vega" }];
    }

    renderPlot() {
        const el = this.container.current;
        if (this.embedded && this.embedded.finalize) {
            this.embedded.finalize();
        }
        const [spec, options] = this.embedArguments();
        return vegaEmbed(el, spec, options).then((embedded) => {
            this.embedded = embedded;
            this.baseWidth = embedded.view.width();
            this.baseHeight = embedded.view.height();
            if (this.fit) {
                return this.fitToWindow().then(() => this.addMenuControls(el));
            }
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

    applyScaleAsync() {
        return this.embedded.view
            .width(this.baseWidth * this.scale / 100)
            .height(this.baseHeight * this.scale / 100)
            .runAsync();
    }

    // scale the plot up or down so that it uses as much of the available space as
    // possible, keeping its aspect ratio
    async fitToWindow() {
        const el = this.container.current;
        const view = this.embedded.view;
        const area = el.parentElement.parentElement;  // the content area of the report
        const drawing = el.querySelector("canvas, svg");
        if (!drawing) {
            return;
        }
        // measure with a tiny plot, because a large plot would itself push the layout
        // and change the available space (no repaint happens in between)
        await view.width(10).height(10).runAsync();
        const padding = 32;
        const availableWidth = area.clientWidth - padding;
        const availableHeight = area.clientHeight - padding;
        // axes, legend and padding around the plot area keep their size when scaling
        const rect = drawing.getBoundingClientRect();
        const extraWidth = rect.width - 10;
        const extraHeight = rect.height - 10;
        const factor = Math.max(
            0.1,
            Math.min(
                (availableWidth - extraWidth) / this.baseWidth,
                (availableHeight - extraHeight) / this.baseHeight
            )
        );
        this.scale = factor * 100;
        await this.applyScaleAsync();
        if (this.sizeSlider) {
            this.sizeSlider.value = Math.round(this.scale / 5) * 5;
            this.sizeValue.textContent = Math.round(this.scale) + "%";
        }
    }

    // the view to save: with "fit to window", the original size
    exportView() {
        if (this.fit) {
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
        menu.append(this.fitRow());
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
        return [row, checkbox];
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
        slider.value = Math.round(this.scale / 5) * 5;
        const value = document.createElement("span");
        value.textContent = Math.round(this.scale) + "%";
        this.sizeSlider = slider;
        this.sizeValue = value;

        // while dragging, only update the label
        slider.addEventListener("input", () => {
            value.textContent = slider.value + "%";
        });
        // re-render once, when the mouse button is released; setting the size by hand
        // switches "fit to window" off
        slider.addEventListener("change", () => {
            this.fit = false;
            this.fitCheckbox.checked = false;
            this.scale = Number(slider.value);
            this.applyScale();
        });

        row.append(label, slider, value);
        return row;
    }

    fitRow() {
        const [row, checkbox] = this.checkboxRow("Fit to window", this.fit, (checked) => {
            this.fit = checked;
            if (checked) {
                this.fitToWindow().catch(this.showError);
            } else {
                // back to the original size
                this.scale = 100;
                this.applyScale();
                this.sizeSlider.value = 100;
                this.sizeValue.textContent = "100%";
            }
        });
        this.fitCheckbox = checkbox;
        return row;
    }

    legendRow() {
        const [row] = this.checkboxRow("Legend", this.showLegend, (checked) => {
            this.showLegend = checked;
            this.renderPlot().catch(this.showError);
        });
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