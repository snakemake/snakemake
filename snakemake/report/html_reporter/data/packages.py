import json
from snakemake.exceptions import WorkflowError
from snakemake.report.html_reporter.common import get_resource_as_string
import snakemake


def get_packages():
    try:
        import pygments
    except ImportError:
        raise WorkflowError(
            "Python package pygments must be installed to create reports."
        )

    return Packages(
        {
            "snakemake": Package(
                version=snakemake.__version__.split("+")[0],
                license_url="assets/snakemake/snakemake/main/LICENSE.md",
            ),
            "pygments": Package(
                version=pygments.__version__,
                license_url="assets/pygments/pygments/master/LICENSE",
            ),
            "tailwindcss": Package(
                version="3.0",
                license_url="assets/tailwindlabs/tailwindcss/3.0/LICENSE",
                url="assets/tailwindlabs/tailwindcss/3.0/3.0.23?plugins=forms@0.4.0,typography@0.5.2"
            ),
            "react": Package(
                version="18.2",
                license_url="assets/facebook/react/18.2/LICENSE",
                main="assets/facebook/react/18.2/react.production.min.js",
                dom="assets/facebook/react/18.2/react-dom.production.min.js",
            ),
            "vega": Package(
                version="5.21",
                url="assets/vega/vega/5.21/vega.js",
                license_url="assets/vega/vega/5.21/LICENSE",
            ),
            "vega-lite": Package(
                version="5.2",
                url="assets/vega/vega-lite/5.2/vega-lite.js",
                license_url="assets/vega/vega-lite/5.2/LICENSE",
            ),
            "vega-embed": Package(
                version="6.20",
                url="assets/vega/vega-embed/6.20/vega-embed.js",
                license_url="assets/vega/vega-embed/6.20/LICENSE",
            ),
            "heroicons": Package(
                version="1.0.6",
                license_url="assets/tailwindlabs/heroicons/master/LICENSE",
            ),
            "prop-types": Package(
                version="15.7.2",
                url="assets/facebook/prop-types/15.7.2/prop-types.min.js",
                license_url="assets/facebook/prop-types/15.7.2/LICENSE",
            ),
        }
    )


class Packages:
    def __init__(self, packages):
        self.packages = packages

    def __getitem__(self, package):
        return self.packages[package]

    def get_json(self):
        return json.dumps(
            {name: package.get_record() for name, package in self.packages.items()}
        )


class Package:
    def __init__(self, version=None, license_url=None, url=None, **urls):
        self.version = version
        self.license = get_resource_as_string(license_url)
        if url is not None:
            self.url = url
        else:
            self.urls = urls

    def get_record(self):
        return {"version": self.version, "license": self.license}
