'use strict';

class Toggle extends React.Component {
    constructor(props) {
        super(props);
        this.state = { value: this.props.defaultValue };
        this.selectValue = this.selectValue.bind(this);

    }

    selectValue(selected) {
        // swap the value
        this.setState({ value: selected });
        this.props.callback(selected);
    }

    render() {
        let props = this.props;
        function valueId(idx) {
            return `toggle-${props.label}-${props.values[idx]}`;
        }

        const radioClasses = "appearance-none hidden";
        const labelClasses = "flex gap-2 items-center transition-all inline-block p-1 has-[:checked]:bg-emerald-500 has-[:checked]:text-slate-800 hover:bg-slate-500 bg-slate-600"

        let selectValue = this.selectValue;

        return e(
            "fieldset",
            {className: "flex gap-2 items-center shrink-0"},
            e(
                "span",
                {className: "uppercase font-bold"},
                `${props.label}`
            ),
            e(
                "span",
                {className: "flex gap-0 shrink-0"},
                e(
                    "label",
                    {for: valueId(0), className: `${labelClasses} rounded-l`},
                    e(
                        "input",
                        {
                            className: radioClasses,
                            type: "radio",
                            checked: this.state.value == props.values[0],
                            id: valueId(0),
                            name: props.label, 
                            value: props.values[0],
                            onClick: () => selectValue(props.values[0])
                        },
                    ),
                    props.values[0]
                ),
                e(
                    "label",
                    {for: valueId(1), className: `${labelClasses} rounded-r`},
                    e(
                        "input",
                        {
                            className: radioClasses,
                            type: "radio",
                            checked: this.state.value == props.values[1], 
                            id: valueId(1),
                            name: props.label, 
                            value: props.values[1], 
                            onClick: () => selectValue(props.values[1])
                        },
                    ),
                    props.values[1]
                )
            )
        )
    }

    // render2() {
    //     return e(
    //         "label",
    //         { for: this.props.label, className: "p-1 relative inline-flex cursor-pointer items-center"},
    //         e(
    //             "span",
    //             {},
    //             this.props.label
    //         ),
    //         e(
    //             "div",
    //             {className: "relative inline-flex cursor-pointer items-center"},
    //             e(
    //                 "input",
    //                 { id: this.props.label, type: "checkbox", className: "peer sr-only", checked: this.state.checked, onChange: this.handleToggle }
    //             ),
    //             e(
    //                 "div",
    //                 { className: "peer flex h-8 items-center gap-4 rounded-full bg-emerald-500 px-3 after:absolute after:left-1 after: after:h-6 after:w-8 after:rounded-full after:bg-white/40 after:transition-all after:content-[''] peer-checked:bg-stone-600 peer-checked:after:translate-x-full peer-focus:outline-none text-sm text-white" },
    //                 e(
    //                     "span",
    //                     {},
    //                     this.props.values[0]
    //                 ),
    //                 e(
    //                     "span",
    //                     {},
    //                     this.props.values[1]
    //                 )
    //             )
    //         )
    //     )
    // }
}