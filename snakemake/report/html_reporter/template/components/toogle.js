'use strict';
import React, { useState } from 'react'

class Toggle extends React.Component {
    constructor(props) {
        super(props);
        this.state = { checked: false };
        this.handleToggle = this.handleToggle.bind(this);
    }

    handleToggle() {
        this.state.checked = !this.state.checked;
    }

    render() {
        return e(
            "label",
            { for: "toggle", className: "flex items-center cursor-pointer select-none text-dark dark:text-white"},
            e(
                "div",
                { className: "relative" },
                e(
                    "input",
                    { id: "toggle", type: "checkbox", className: "peer sr-only", onChange: this.handleToggle }
                ),
                e(
                    "div",
                    { className: "'block h-8 w-14 rounded-full bg-emerald-500" }
                ),
                e(
                    "div",
                    { className: "dot absolute left-1 top-1 h-6 w-6 rounded-full bg-white transition" }
                )
            )
        )
    }
}