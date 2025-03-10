'use strict';

class AbstractMenu extends React.Component {
    buttonProps = { className: "transition-all block hover:text-emerald-500 rounded hover:bg-slate-800 p-1 flex items-center gap-2" };
    iconProps = { className: "text-emerald-500" };

    getMenuItem(label, iconName, onClick) {
        return e(
            "li",
            { key: label },
            e(
                "a",
                { href: "#", onClick: onClick, ...this.buttonProps },
                e(Icon, { iconName: iconName, ...this.iconProps }),
                e(
                    "span",
                    {},
                    label
                )
            )
        );
    }
}