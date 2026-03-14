'use strict';

class ReportInfo extends AbstractMenu {
    render() {
        return e(
            "ul",
            {},
            e(
                ListHeading,
                { text: "Embedded packages" },
            ),
            this.getPackages()
        );
    }

    getPackages() {
        let _this = this;
        return Object.entries(packages).map(function ([name, record]) {
            return _this.getMenuItem(`${name} ${record.version}`, "cube", () => _this.props.app.showLicense(name))
        }
        );
    }
}