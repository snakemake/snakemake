'use strict';

class VegaPlot extends React.Component {
    constructor(props) {
        super(props);
        this.container = React.createRef();
        this.embedded = null;
        this.spec = null;
        this.isVegaLite = false;
        this.hasLegend = false;
        this.sizeSignals = [];
        this.fixedWidth = null;
        this.fixedHeight = null;
        this.scale = 100;
        this.fitScale = 100;  // the scale at which the plot fits the window
        this.fontScale = 100;
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
        if (!this.embedded || this.resizeFrame !== null) {
            return;
        }
        this.resizeFrame = requestAnimationFrame(() => {
            this.resizeFrame = null;
            this.updateScale().catch(this.showError);
        });
    }

    // plain Vega version of a spec (a copy)
    compiledSpec(spec) {
        return this.isVegaLite
            ? vegaLite.compile(spec, { config: this.fontConfig() }).spec
            : JSON.parse(JSON.stringify(spec));
    }

    // font sizes of the plot (the Vega-Lite defaults, multiplied by the font scale)
    fontConfig() {
        const f = this.fontScale / 100;
        return {
            axis: { labelFontSize: 10 * f, titleFontSize: 11 * f },
            legend: { labelFontSize: 10 * f, titleFontSize: 11 * f },
            title: { fontSize: 13 * f, subtitleFontSize: 12 * f },
            header: { labelFontSize: 10 * f, titleFontSize: 10 * f },
            text: { fontSize: 11 * f }
        };
    }

    // spec and options for the current settings (legend shown or hidden, font size)
    embedArguments() {
        const options = {
            actions: { export: false, source: true, compiled: true, editor: true },
            forceActionsMenu: true,
            config: this.fontConfig()
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
        const details = el.querySelector("details");
        const menuWasOpen = details !== null && details.hasAttribute("open");
        const [spec, options] = this.embedArguments();
        return vegaEmbed(el, spec, options).then((embedded) => {
            this.embedded = embedded;
            this.collectSizes(embedded.vgSpec);
            const finish = () => {
                this.addMenuControls(el);
                if (menuWasOpen) {
                    el.querySelector("details").setAttribute("open", "");
                }
            };
            return this.updateScale().then(finish);
        });
    }

    // Everything that sets the size of the plot area in the compiled Vega spec:
    // plain width and height, and the sizes of discrete axes ("step") and of
    // facets or concatenated plots (e.g. x_step, child_width, concat_0_width)
    collectSizes(vgSpec) {
        this.sizeSignals = (vgSpec.signals || [])
            .filter((signal) => typeof signal.value === "number" && /(width|height|step)$/i.test(signal.name))
            .map((signal) => ({ name: signal.name, value: signal.value }));
        this.fixedWidth = typeof vgSpec.width === "number" ? vgSpec.width : null;
        this.fixedHeight = typeof vgSpec.height === "number" ? vgSpec.height : null;
    }

    // set the size of the plot area to factor times the original size
    setScale(factor) {
        const view = this.embedded.view;
        for (const signal of this.sizeSignals) {
            view.signal(signal.name, signal.value * factor);
        }
        if (this.fixedWidth !== null) {
            view.width(this.fixedWidth * factor);
        }
        if (this.fixedHeight !== null) {
            view.height(this.fixedHeight * factor);
        }
        // faceted and concatenated plots only change their size after a resize
        return view.resize().runAsync();
    }

    // Find the scale at which the plot uses as much of the available space as
    // possible, keeping its aspect ratio (the plot is left at that scale)
    async computeFitFactor() {
        const el = this.container.current;
        const area = el.parentElement.parentElement;  // the content area of the report
        const drawing = el.querySelector("canvas, svg");
        if (!drawing) {
            return 1;
        }
        const measure = () => {
            const rect = drawing.getBoundingClientRect();
            return [rect.width, rect.height];
        };
        // Axes, legend and padding keep their size when the plot area is scaled, so
        // the total size is roughly (fixed part) + factor * (scaled part). Measure it
        // at the original size and at a small size (this also keeps a large plot from
        // pushing the layout while the available space is measured), then correct the
        // estimate a few times, because e.g. axis labels can change with the size.
        const points = [];
        for (const factor of [1, 0.2]) {
            await this.setScale(factor);
            points.push({ factor, size: measure() });
        }
        const padding = 32;
        const available = [area.clientWidth - padding, area.clientHeight - padding];
        let factor = 1;
        for (let i = 0; i < 4; i++) {
            const [a, b] = points.slice(-2);
            const estimates = [0, 1].map((dim) => {
                const slope = (a.size[dim] - b.size[dim]) / (a.factor - b.factor);
                return slope > 0.5 ? a.factor + (available[dim] - a.size[dim]) / slope : Infinity;
            });
            factor = Math.max(0.1, Math.min(...estimates));
            if (!isFinite(factor)) {
                factor = 1;
                break;
            }
            await this.setScale(factor);
            const size = measure();
            points.push({ factor, size });
            const fits = size[0] <= available[0] + 1 && size[1] <= available[1] + 1;
            const full = size[0] >= available[0] - 4 || size[1] >= available[1] - 4;
            if (fits && full) {
                break;
            }
        }
        return factor;
    }

    // Update the scale at which the plot fits the window, which is the maximum of the
    // size slider. With "fit to window" the plot is set to that scale, otherwise it
    // keeps its scale (but never gets larger than the fitting one).
    async updateScale() {
        this.fitScale = (await this.computeFitFactor()) * 100;
        if (this.fit) {
            this.scale = this.fitScale;
        } else {
            this.scale = Math.min(this.scale, this.fitScale);
            await this.setScale(this.scale / 100);
        }
        this.updateSizeSlider();
    }

    updateSizeSlider() {
        if (!this.sizeSlider) {
            return;
        }
        this.sizeSlider.min = Math.min(25, this.fitScale);
        this.sizeSlider.max = this.fitScale;
        this.sizeSlider.value = this.scale;
        this.sizeValue.textContent = Math.round(this.scale) + "%";
    }

    // the view to save: with "fit to window", the original size
    exportView() {
        if (this.fit) {
            const vgSpec = this.compiledSpec(this.spec);
            if (!this.showLegend) {
                delete vgSpec.legends;
            }
            const view = new vega.View(vega.parse(vgSpec, this.fontConfig()), { renderer: "none" });
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
        menu.append(this.fontRow());
        menu.append(this.fitRow());
        if (this.hasLegend) {
            menu.append(this.legendRow());
        }
        menu.append(this.restoreItem());
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
        // the maximum is the scale at which the plot fits the window
        slider.min = Math.min(25, this.fitScale);
        slider.max = this.fitScale;
        slider.step = "any";
        slider.value = this.scale;
        const value = document.createElement("span");
        value.textContent = Math.round(this.scale) + "%";
        this.sizeSlider = slider;
        this.sizeValue = value;

        // while dragging, only update the label
        slider.addEventListener("input", () => {
            value.textContent = Math.round(slider.value) + "%";
        });
        // re-render once, when the mouse button is released; setting the size by hand
        // switches "fit to window" off
        slider.addEventListener("change", () => {
            this.fit = false;
            this.fitCheckbox.checked = false;
            this.scale = Number(slider.value);
            this.setScale(this.scale / 100);
        });

        row.append(label, slider, value);
        return row;
    }

    fontRow() {
        const row = this.menuRow();
        const label = document.createElement("span");
        label.textContent = "Font";
        const slider = document.createElement("input");
        slider.type = "range";
        slider.min = 50;
        slider.max = 300;
        slider.step = 10;
        slider.value = this.fontScale;
        const value = document.createElement("span");
        value.textContent = this.fontScale + "%";
        slider.addEventListener("input", () => {
            value.textContent = slider.value + "%";
        });
        // the fonts are changed once, when the mouse button is released
        slider.addEventListener("change", () => {
            this.fontScale = Number(slider.value);
            this.renderPlot().catch(this.showError);
        });
        row.append(label, slider, value);
        return row;
    }

    fitRow() {
        const [row, checkbox] = this.checkboxRow("Fit to window", this.fit, (checked) => {
            this.fit = checked;
            if (checked) {
                this.scale = this.fitScale;
                this.setScale(this.scale / 100);
                this.updateSizeSlider();
            }
        });
        this.fitCheckbox = checkbox;
        return row;
    }
    // back to the defaults: fit to window, normal fonts, legend shown
    restoreItem() {
        const item = document.createElement("a");
        item.textContent = "Restore defaults";
        item.addEventListener("click", () => {
            this.fit = true;
            this.fontScale = 100;
            this.showLegend = true;
            this.renderPlot().catch(this.showError);
        });
        return item;
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