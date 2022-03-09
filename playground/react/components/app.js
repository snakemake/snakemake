'use strict';



class App extends React.Component {
    constructor(props) {
        super(props);
        this.state = { navbarMode: "menu", content: "rulegraph", ruleinfo: undefined, category: undefined, subcategory: undefined, searchTerm: undefined, resultPath: undefined };
        this.setView = this.setView.bind(this);
    }

    render() {
        return [
            e(Navbar, { app: this }),
            e(ContentDisplay, { show: this.state.content })
        ];
    }

    setView(view) {
        event.preventDefault();
        this.setState({
            navbarMode: view.mode || this.state.navbarMode,
            content: view.content || this.state.content,
            ruleinfo: view.ruleinfo,
            category: view.category,
            subcategory: view.subcategory,
            searchTerm: view.searchTerm,
            resultPath: view.resultPath
        })
    }
}

ReactDOM.render(e(App), document.querySelector('#app'));