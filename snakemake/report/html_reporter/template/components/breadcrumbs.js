'use strict';

class Breadcrumbs extends React.Component {
    constructor(props) {
        super(props);
        this.state = { searchTerm: "" };
    }

    render() {
        return e(
            "nav",
            { className: "text-white whitespace-nowrap align-middle text-xs m-2 p-2 rounded bg-slate-800 min-w-fit" },
            e(
                "ol",
                { className: "list-reset flex items-center gap-1" },
                this.renderEntries(),
                this.renderSearch()
            )
        )
    }

    renderEntries() {
        let entries = this.props.entries.filter(function (entry) {
            return entry !== undefined
        });

        return entries.map(function (entry, index) {
            const isLast = index == entries.length - 1;

            let content = entry.name;
            if (entry.icon !== undefined) {
                content = e(
                    Icon,
                    { iconName: entry.icon }
                )
            }


            let link = content;
            if (entry.func !== undefined) {
                link = e(
                    "a",
                    { className: "hover:text-emerald-600", href: "#", onClick: entry.func },
                    content
                );
            }

            let props = {};
            if (isLast) {
                props = { className: "grow" }
            }

            let item = e(
                "li",
                { key: entry.name, ...props },
                link
            );
            if (!isLast) {
                item = [
                    item,
                    e(
                        "li",
                        { key: `sep-${index}` },
                        e(
                            Icon,
                            { iconName: "chevron-right", className: "text-emerald-500" }
                        )
                    )
                ];
            }

            return item;
        });
    }

    renderSearch() {
        let _this = this;
        return e(
            "li",
            { className: "flex-0" },
            e(
                "input",
                {
                    type: "text",
                    placeholder: "Search...",
                    title: "Search all results. Prefix with 're:' to perform regexp search.",
                    size: 10,
                    className: "border-0 bg-transparent text-white h-3 w-fit text-xs text-right form-control rounded",
                    onKeyPress: function (event) {
                        if (event.charCode == 13) {
                            event.preventDefault();
                            _this.showSearchResults();
                        }
                    },
                    onChange: function (event) {
                        _this.setState({ searchTerm: event.target.value });
                    }
                },
            )
        )
    }

    showSearchResults() {
        this.props.setView({ navbarMode: "searchresults", searchTerm: this.state.searchTerm })
    }
}

Breadcrumbs.propTypes = {
    entries: PropTypes.arrayOf(PropTypes.object).isRequired,
    setView: PropTypes.func.isRequired,
};