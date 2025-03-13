'use strict';

class Category extends AbstractMenu {
    render() {
        return e(
            "ul",
            {},
            this.getSubcategoryMenuitems()
        )
    }

    showSubcategory(subcategory) {
        this.props.setView({ navbarMode: "subcategory", category: this.props.category, subcategory: subcategory })
    }

    getSubcategoryMenuitems() {
        let _this = this
        let items = Object.keys(categories[this.props.category]).sort(
            (a, b) => a.localeCompare(b)
        ).map(function (subcategory) {
            return _this.getMenuItem(subcategory, "folder", () => _this.showSubcategory(subcategory));
        });
        return items;
    }
}