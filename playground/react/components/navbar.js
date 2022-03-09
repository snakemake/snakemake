'use strict';


class Navbar extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        let translateNavbar = "";
        let translateShowButton = "-translate-x-full"
        if (this.props.app.state.hideNavbar) {
            translateNavbar = "-translate-x-full";
            translateShowButton = ""
        }
        return [
            e(
                "div",
                { className: `fixed z-50 p-3 transition-translate ${translateShowButton}` },
                this.getShowButton()
            ),
            e(
                "nav",
                { className: `fixed z-50 transition-all ${translateNavbar} ${this.getWidth()} min-w-fit text-white text-sm bg-slate-900/70 backdrop-blur-sm h-screen overflow-y-auto` },
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
                        { className: "grow mr-5" },
                        "Report"
                    ),
                    this.getHideButton()
                ),
                this.renderBreadcrumbs(),
                e(
                    "div",
                    { className: "p-3" },
                    this.renderContent()
                )
            )
        ];
    }

    getHideButton() {
        let setView = this.props.app.setView;
        return e(
            "a",
            {
                type: "button",
                className: "bg-transparent hover:text-emerald-500",
                href: "#",
                onClick: () => setView({ hideNavbar: true })
            },
            e(
                Icon,
                { iconName: "x" }
            )
        )
    }

    getShowButton() {
        let setView = this.props.app.setView;
        return e(
            "a",
            {
                type: "button",
                href: "#",
                className: "bg-transparent hover:text-emerald-500",
                onClick: () => setView({ hideNavbar: false })
            },
            e(
                Icon,
                { iconName: "menu" }
            )
        )
    }

    renderContent() {
        let setView = this.props.app.setView;
        switch (this.props.app.state.navbarMode) {
            case "menu":
                return e(Menu, { setView: setView, app: this.props.app });
            case "category":
                return e(Category, { setView: setView, category: this.props.app.state.category });
            case "subcategory":
                return e(Subcategory, { setView: setView, category: this.props.app.state.category, subcategory: this.props.app.state.subcategory });
            case "searchresults":
                return e(SearchResults, { setView: setView, searchTerm: this.props.app.state.searchTerm });
            case "resultinfo":
                return e(ResultInfo, { resultPath: this.props.app.state.resultPath, setView: setView });
            case "ruleinfo":
                return e(RuleInfo, { rule: this.props.app.state.ruleinfo });
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
                    { entries: [this.getMenuBreadcrumb(), this.getResultBreadcrumb(), this.getCategoryBreadcrumb()], setView: setView }
                );
            case "subcategory":
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
            case "ruleinfo":
                if (this.props.app.state.resultPath) {
                    return e(
                        Breadcrumbs,
                        { entries: [this.getMenuBreadcrumb(), this.getResultBreadcrumb(), this.getCategoryBreadcrumb(), this.getSubcategoryBreadcrumb(), this.getSearchResultsBreadcrumb(), this.getResultinfoBreadcrumb(), this.getRuleBreadcrumb(), this.getRuleinfoBreadcrumb()], setView: setView }
                    )
                } else {
                    return e(
                        Breadcrumbs,
                        { entries: [this.getMenuBreadcrumb(), this.getRuleBreadcrumb(), this.getRuleinfoBreadcrumb()], setView: setView }
                    )
                }
        }
    }

    getMenuBreadcrumb() {
        let setView = this.props.app.setView;
        return { name: "menu", icon: "home", func: function () { setView({ navbarMode: "menu", category: undefined, subcategory: undefined }) } };
    }

    getRuleBreadcrumb() {
        return { name: "Rule", func: undefined }
    }

    getRuleinfoBreadcrumb() {
        return { name: this.props.app.state.ruleinfo, func: undefined }
    }

    getCategoryBreadcrumb() {
        let category = this.props.app.state.category;
        if (category === undefined) {
            return undefined;
        }
        let subcategory = undefined;
        let mode = "category";
        if (isSingleSubcategory(category)) {
            subcategory = this.props.app.state.subcategory;
            mode = "subcategory";
        }
        let _this = this;

        let name = this.props.app.state.category;
        if (isSingleDefaultCategory()) {
            name = "Results";
        }
        let setView = this.props.app.setView;
        return { name: name, func: function () { setView({ navbarMode: mode, category: category, subcategory: subcategory }) } };
    }

    getSubcategoryBreadcrumb() {
        let subcategory = this.props.app.state.subcategory;
        let category = this.props.app.state.category;
        if (subcategory === undefined || isSingleSubcategory(category)) {
            return undefined;
        }
        let setView = this.props.app.setView;
        return { name: this.props.app.state.subcategory, func: function () { setView({ navbarMode: "subcategory", category: category, subcategory: subcategory }) } };
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
                setView({ navbarMode: "searchresults", searchTerm: searchTerm })
            }
        };
    }

    getResultinfoBreadcrumb() {
        let setView = this.props.app.setView;
        return {
            name: "resultinfo", icon: "eye", func: function () { setView({ navbarMode: "resultinfo" }) }
        };
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