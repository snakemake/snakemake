'use strict';

class Menu extends AbstractMenu {
    constructor(props) {
        super(props);
        this.showWorkflow = this.showWorkflow.bind(this)
        this.showStatistics = this.showStatistics.bind(this)
    }

    render() {
        return e(
            "ul",
            {},
            this.getHeading(),
            this.getMenuItem("Workflow", "share", this.showWorkflow),
            this.getMenuItem("Statistics", "chart", this.showStatistics),
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

    showCategory(category) {
        let subcategory = undefined;
        let mode = "category";
        if (isSingleSubcategory(category)) {
            subcategory = Object.keys(categories[category])[0];
            mode = "subcategory";
        }
        this.props.setView({ navbarMode: mode, category: category, subcategory: subcategory })
    }

    getCategoryMenumitems() {
        if (isSingleCategory()) {
            let category = Object.keys(categories)[0];
            return this.getMenuItem("Results", "folder", () => this.showCategory(category));
        } else if (isNoResults()) {
            return [];
        } else {
            let items = [e(
                "li",
                { key: "Results", className: "uppercase font-bold p-1 mt-2" },
                "Results"
            )];

            let _this = this
            items.push(...Object.keys(categories).map(function (category) {
                return _this.getMenuItem(category, "folder", () => _this.showCategory(category));
            }));

            return items;
        }
    }
}