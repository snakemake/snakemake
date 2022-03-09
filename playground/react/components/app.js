'use strict';



class App extends React.Component {
    constructor(props) {
        super(props);
        this.state = { content: "rulegraph", ruleinfo: undefined };
    }

    render() {
        return [
            e(Navbar, { app: this, ruleinfo: this.state.ruleinfo }),
            e(ContentDisplay, { show: this.state.content })
        ];
    }

    setView(view) {
        this.setState({ content: view.content, ruleinfo: view.ruleinfo })
    }
}

ReactDOM.render(e(App), document.querySelector('#app'));