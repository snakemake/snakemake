'use strict';

class Toggle extends React.Component {
    constructor(props) {
        super(props);
        this.state = { value: this.props.values[0], checked: false };
        this.handleToggle = this.handleToggle.bind(this);
    }

    handleToggle() {
        // swap the value
        let value = this.props.values.find((value) => value !== this.state.value);
        this.setState({ value: value, checked: !this.state.checked });
        this.props.callback(value);
    }

    render() {
        return e(
            "label",
            { for: this.props.label, className: "p-1 relative inline-flex cursor-pointer items-center"},
            e(
                "span",
                {},
                this.props.label
            ),
            e(
                "div",
                {className: "relative inline-flex cursor-pointer items-center"},
                e(
                    "input",
                    { id: this.props.label, type: "checkbox", className: "peer sr-only", checked: this.state.checked, onChange: this.handleToggle }
                ),
                e(
                    "div",
                    { className: "peer flex h-8 items-center gap-4 rounded-full bg-emerald-500 px-3 after:absolute after:left-1 after: after:h-6 after:w-8 after:rounded-full after:bg-white/40 after:transition-all after:content-[''] peer-checked:bg-stone-600 peer-checked:after:translate-x-full peer-focus:outline-none text-sm text-white" },
                    e(
                        "span",
                        {},
                        this.props.values[0]
                    ),
                    e(
                        "span",
                        {},
                        this.props.values[1]
                    )
                )
            )
        )
    }
}