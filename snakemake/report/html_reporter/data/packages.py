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
                license_url="https://raw.githubusercontent.com/snakemake/snakemake/main/LICENSE.md",
            ),
            "pygments": Package(
                version=pygments.__version__,
                license_url="https://raw.githubusercontent.com/pygments/pygments/master/LICENSE",
            ),
            "tailwindcss": Package(
                version="3.0",
                license_url="https://raw.githubusercontent.com/tailwindlabs/tailwindcss/master/LICENSE",
                url="https://cdn.tailwindcss.com/3.0.23?plugins=forms@0.4.0,typography@0.5.2",
            ),
            "react": Package(
                version="18",
                license_url="https://raw.githubusercontent.com/facebook/react/main/LICENSE",
                main="https://cdnjs.cloudflare.com/ajax/libs/react/18.2.0/umd/react.production.min.js",
                dom="https://cdnjs.cloudflare.com/ajax/libs/react-dom/18.2.0/umd/react-dom.production.min.js",
            ),
            "vega": Package(
                version="5.21",
                url="https://cdnjs.cloudflare.com/ajax/libs/vega/5.21.0/vega.js",
                license_url="https://raw.githubusercontent.com/vega/vega/main/LICENSE",
            ),
            "vega-lite": Package(
                version="5.2",
                url="https://cdnjs.cloudflare.com/ajax/libs/vega-lite/5.2.0/vega-lite.js",
                license_url="https://raw.githubusercontent.com/vega/vega-lite/next/LICENSE",
            ),
            "vega-embed": Package(
                version="6.20",
                url="https://cdnjs.cloudflare.com/ajax/libs/vega-embed/6.20.8/vega-embed.js",
                license_url="https://raw.githubusercontent.com/vega/vega-embed/next/LICENSE",
            ),
            "heroicons": Package(
                version="1.0.6",
                license_url="https://raw.githubusercontent.com/tailwindlabs/heroicons/master/LICENSE",
            ),
            "prop-types": Package(
                version="15.7.2",
                url="https://cdnjs.cloudflare.com/ajax/libs/prop-types/15.7.2/prop-types.min.js",
                license_url="https://raw.githubusercontent.com/facebook/prop-types/main/LICENSE",
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
