'use strict';


let app;


class App extends React.Component {
    constructor(props) {
        super(props);
        this.state = { hideNavbar: false, navbarMode: "menu", content: "rulegraph", ruleinfo: undefined, category: undefined, subcategory: undefined, searchTerm: undefined, resultPath: undefined, contentPath: undefined };
        this.setView = this.setView.bind(this);
        this.showCategory = this.showCategory.bind(this);
        this.showResultInfo = this.showResultInfo.bind(this);
        this.showReportInfo = this.showReportInfo.bind(this);
        // store in global variable
        app = this;
    }

    render() {
        return [
            e(
                "div",
                { class: "flex flex-row w-screen h-screen" },
                e(Navbar, { key: "navbar", app: this }),
                e(ContentDisplay, { key: "content", app: this })
            )
        ];
    }

    setView(view) {
        this.setState({
            hideNavbar: view.hideNavbar || this.hideNavbar,
            navbarMode: view.navbarMode || this.state.navbarMode,
            content: view.content || this.state.content,
            ruleinfo: view.ruleinfo || this.state.ruleinfo,
            category: view.category || this.state.category,
            subcategory: view.subcategory || this.state.subcategory,
            searchTerm: view.searchTerm || this.state.searchTerm,
            resultPath: view.resultPath || this.state.resultPath,
            contentPath: view.contentPath || this.state.contentPath,
            contentText: view.contentText || this.state.contentText,
        })
    }

    showCategory(category) {
        let subcategory;
        let mode = "category";
        if (isSingleSubcategory(category)) {
            subcategory = Object.keys(categories[category])[0];
            mode = "subcategory";
        }
        this.setView({ navbarMode: mode, category: category, subcategory: subcategory })
    }

    showReportInfo() {
        this.setView({ navbarMode: "reportinfo" });
    }

    showResultInfo(resultPath) {
        this.setView({ navbarMode: "resultinfo", resultPath: resultPath });
    }

    showLicense(package_name) {
        this.setView({
            content: "text",
            contentText: packages[package_name].license
        });
    }
}

const root = ReactDOM.createRoot(document.querySelector('#app'));
root.render(e(App));