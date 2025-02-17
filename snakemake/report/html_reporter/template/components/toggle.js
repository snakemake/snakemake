'use strict';

class Toggle extends React.Component {
    constructor(props) {
        super(props);
        this.state = { value: this.props.values[0] };
        this.selectValue = this.selectValue.bind(this);

    }

    selectValue(selected) {
        // swap the value
        console.log(selected);
        this.setState({ value: selected });
        this.props.callback(selected);
    }

    render() {
        let props = this.props;
        function valueId(idx) {
            return `toggle-${props.label}-${props.values[idx]}`;
        }

        const radioClasses = "focus:shadow-none outline-none focus:ring-0 focus:outline-none appearance-none w-4 h-4 checked:bg-emerald-500 focus:shadow-none focus:ring-emerald-500 focus:bg-slate-500 focus:border-emerald-500 focus:border-2 outline-none bg-clip-padding border-2 border-emerald-500 rounded-full shrink-0 bg-slate-500";

        let selectValue = this.selectValue;

        return e(
            "fieldset",
            {className: "flex gap-2"},
            e(
                "span",
                {},
                `${props.label}: `
            ),
            e(
                "span",
                {className: "flex gap-2 items-center"},
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
                e(
                    "label",
                    {for: valueId(0)},
                    props.values[0]
                )
            ),
            e(
                "span",
                {className: "flex gap-2 items-center"},
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
                e(
                    "label",
                    {for: valueId(1)},
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