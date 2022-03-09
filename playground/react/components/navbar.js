'use strict';


class Navbar extends React.Component {
    constructor(props) {
        super(props);
        this.state = { mode: "menu", category: undefined, subcategory: undefined, searchTerm: undefined, resultPath: undefined };
        this.setView = this.setView.bind(this);
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
        let setView = this.setView;
        switch (this.state.mode) {
            case "menu":
                return e(Menu, { setView: setView, app: this.props.app });
            case "category":
                if (this.state.subcategory !== undefined) {
                    return e(Subcategory, { setView: setView, category: this.state.category, subcategory: this.state.subcategory });
                } else {
                    return e(Category, { setView: setView, category: this.state.category });
                }
            case "searchresults":
                return e(SearchResults, { setView: setView, searchTerm: this.state.searchTerm });
            case "resultinfo":
                return e(ResultInfo, { resultPath: this.state.resultPath });
        }
    }

    renderBreadcrumbs() {
        let setView = this.setView;
        switch (this.state.mode) {
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
        let _this = this;
        return { name: "menu", icon: "home", func: function () { _this.setView({ mode: "menu", category: undefined, subcategory: undefined }) } };
    }

    getCategoryBreadcrumb() {
        if (this.state.category === undefined) {
            return undefined;
        }
        let subcategory = undefined;
        if (isSingleSubcategory(this.state.category)) {
            subcategory = this.state.subcategory;
        }
        let _this = this;

        let name = this.state.category;
        if (isSingleDefaultCategory()) {
            name = "Results";
        }
        return { name: name, func: function () { _this.setView({ mode: "category", category: _this.state.category, subcategory: subcategory }) } };
    }

    getSubcategoryBreadcrumb() {
        if (this.state.subcategory === undefined || isSingleSubcategory(this.state.category)) {
            return undefined;
        }
        let _this = this;
        return { name: this.state.subcategory, func: function () { _this.setView({ mode: "category", category: _this.state.category, subcategory: _this.state.subcategory }) } };
    }

    getResultBreadcrumb() {
        if (isSingleDefaultCategory()) {
            return undefined;
        }
        return { name: "Results", func: undefined };
    }

    getSearchResultsBreadcrumb() {
        let searchTerm = this.state.searchTerm;
        if (searchTerm === undefined) {
            return undefined;
        }
        let setView = this.setView;
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
        switch (this.state.mode) {
            case "menu":
                return "w-1/5"
            case "category":
                return "w-1/3"
            case "resultinfo":
                return "w-3/4"
        }
    }

    setView(view) {
        event.preventDefault();
        this.setState({ mode: view.mode, category: view.category, subcategory: view.subcategory, searchTerm: view.searchTerm, resultPath: view.resultPath })
    }
}