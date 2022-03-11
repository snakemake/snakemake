'use strict';



class App extends React.Component {
    constructor(props) {
        super(props);
        this.state = { hideNavbar: false, navbarMode: "menu", content: "rulegraph", ruleinfo: undefined, category: undefined, subcategory: undefined, searchTerm: undefined, resultPath: undefined, contentPath: undefined };
        this.setView = this.setView.bind(this);
    }

    render() {
        return [
            e(Navbar, { key: "navbar", app: this }),
            e(ContentDisplay, { key: "content", app: this })
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
        })
    }
}

ReactDOM.render(e(App), document.querySelector('#app'));