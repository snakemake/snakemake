'use strict';

class Menu extends AbstractMenu {
    constructor(props) {
        super(props);
        this.showWorkflow = this.showWorkflow.bind(this);
        this.showStatistics = this.showStatistics.bind(this);
    }

    render() {
        return e(
            "ul",
            {},
            this.getHeading(),
            this.getMenuItem("Workflow", "share", this.showWorkflow),
            this.getMenuItem("Statistics", "chart", this.showStatistics),
            this.getMenuItem("About", "information-circle", this.props.app.showReportInfo),
            this.getCategoryMenumitems()
        )
    }

    showWorkflow() {
        this.props.app.setView({ content: "rulegraph" });
    }

    showStatistics() {
        this.props.app.setView({ content: "stats" });
    }

    getHeading() {
        if (isSingleCategory()) {
            return [];
        } else {
            return e(
                ListHeading,
                { text: "General" }
            )
        }
    }

    getCategoryMenumitems() {
        let _this = this;
        let app = this.props.app;
        if (isSingleCategory()) {
            let category = Object.keys(categories)[0];
            return this.getMenuItem("Results", "folder", () => app.showCategory(category));
        } else if (isNoResults()) {
            return [];
        } else {
            let items = [e(
                ListHeading,
                { key: "Results", text: "Results" }
            )];

            items.push(...Object.keys(categories).sort(
                (a, b) => a.localeCompare(b)
            ).map(function (category) {
                return _this.getMenuItem(category, "folder", () => app.showCategory(category));
            }));

            return items;
        }
    }
}