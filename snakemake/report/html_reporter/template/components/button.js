'use strict';

class Button extends React.Component {
    
    render() {
        return this.renderButton(this.props.iconName, this.props);
    }

    renderButton(iconName, props) {
        return e(
            "a",
            { type: "button", className: `transition-all inline-block p-1 text-emerald-500 rounded hover:bg-slate-800`, ...props },
            e(Icon, { iconName: iconName })
        );
    }
}

Button.propTypes = {
    iconName: PropTypes.string.isRequired,
};
