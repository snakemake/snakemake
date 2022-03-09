'use strict';


class Navbar extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        return e(
            "nav",
            { className: `fixed z-50 transition-all ${this.getWidth()} min-w-fit bg-slate-900/70 backdrop-blur-sm text-white text-sm h-screen overflow-y-auto` },
            e(
                "h1",
                { className: "sticky bg-blur bg-white opacity-80 text-slate-700 text-l tracking-wide px-3 py-1 mb-1 flex items-center" },
                e(
                    "img",
                    { src: "logo.svg", className: "h-4" }
                ),
                e(
                    "span",
                    { className: "font-bold mx-1" },
                    "Snakemake"
                ),
                e(
                    "span",
                    {},
                    "Report"
                ),
            ),
            this.renderBreadcrumbs(),
            e(
                "div",
                { className: "p-3" },
                this.renderContent()
            )
        );
    }

    renderContent() {
        let setView = this.props.app.setView;
        switch (this.props.app.state.navbarMode) {
            case "menu":
                return e(Menu, { setView: setView, app: this.props.app });
            case "category":
                if (this.props.app.state.subcategory !== undefined) {
                    return e(Subcategory, { setView: setView, category: this.props.app.state.category, subcategory: this.props.app.state.subcategory });
                } else {
                    return e(Category, { setView: setView, category: this.props.app.state.category });
                }
            case "searchresults":
                return e(SearchResults, { setView: setView, searchTerm: this.props.app.state.searchTerm });
            case "resultinfo":
                return e(ResultInfo, { resultPath: this.props.app.state.resultPath });
        }
    }

    renderBreadcrumbs() {
        let setView = this.props.app.setView;
        switch (this.props.app.state.navbarMode) {
            case "menu":
                return e(
                    Breadcrumbs,
                    { entries: [this.getMenuBreadcrumb()], setView: setView }
                );
            case "category":
                return e(
                    Breadcrumbs,
                    { entries: [this.getMenuBreadcrumb(), this.getResultBreadcrumb(), this.getCategoryBreadcrumb(), this.getSubcategoryBreadcrumb()], setView: setView }
                );
            case "resultinfo":
                return e(
                    Breadcrumbs,
                    { entries: [this.getMenuBreadcrumb(), this.getResultBreadcrumb(), this.getCategoryBreadcrumb(), this.getSubcategoryBreadcrumb(), this.getSearchResultsBreadcrumb(), this.getResultinfoBreadcrumb()], setView: setView }
                );
            case "searchresults":
                return e(
                    Breadcrumbs,
                    { entries: [this.getMenuBreadcrumb(), this.getSearchResultsBreadcrumb()], setView: setView }
                )
        }
    }

    getMenuBreadcrumb() {
        let setView = this.props.app.setView;
        return { name: "menu", icon: "home", func: function () { setView({ mode: "menu", category: undefined, subcategory: undefined }) } };
    }

    getCategoryBreadcrumb() {
        let category = this.props.app.state.category;
        if (category === undefined) {
            return undefined;
        }
        let subcategory = undefined;
        if (isSingleSubcategory(category)) {
            subcategory = this.props.app.state.subcategory;
        }
        let _this = this;

        let name = this.props.app.state.category;
        if (isSingleDefaultCategory()) {
            name = "Results";
        }
        let setView = this.props.app.setView;
        return { name: name, func: function () { setView({ mode: "category", category: category, subcategory: subcategory }) } };
    }

    getSubcategoryBreadcrumb() {
        let subcategory = this.props.app.state.subcategory;
        let category = this.props.app.state.category;
        if (subcategory === undefined || isSingleSubcategory(category)) {
            return undefined;
        }
        let setView = this.props.app.setView;
        return { name: this.props.app.state.subcategory, func: function () { setView({ mode: "category", category: category, subcategory: subcategory }) } };
    }

    getResultBreadcrumb() {
        if (isSingleDefaultCategory()) {
            return undefined;
        }
        return { name: "Results", func: undefined };
    }

    getSearchResultsBreadcrumb() {
        let searchTerm = this.props.app.state.searchTerm;
        if (searchTerm === undefined) {
            return undefined;
        }
        let setView = this.props.app.setView;
        return {
            name: "Search results", func: function () {
                setView({ mode: "searchresults", searchTerm: searchTerm })
            }
        };
    }

    getResultinfoBreadcrumb() {
        return { name: "resultinfo", icon: "eye", func: undefined };
    }

    getWidth() {
        switch (this.props.app.state.navbarMode) {
            case "menu":
                return "w-1/5"
            case "category":
                return "w-1/3"
            case "resultinfo":
                return "w-3/4"
        }
    }
}