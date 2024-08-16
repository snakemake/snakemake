import json
from snakemake.exceptions import WorkflowError
from snakemake.report.html_reporter.common import get_resource_as_string
import snakemake
import sys
import os


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
                license_url=os.path.join(sys.prefix, "share/snakemake/LICENSE.md"),
            ),
            "pygments": Package(
                version=pygments.__version__,
                license_url=os.path.join(sys.prefix, "share/pygments/LICENSE"),
            ),
            "tailwindcss": Package(
                license_url=os.path.join(sys.prefix, "share/tailwindcss/LICENSE"),
                url=os.path.join(sys.prefix, "share/tailwindcss/tailwind.css"),
            ),
            "react": Package(
                license_url=os.path.join(sys.prefix, "share/react/LICENSE"),
                main=os.path.join(sys.prefix, "share/react/react.production.min.js"),
                dom=os.path.join(sys.prefix, "share/react/react-dom.production.min.js"),
            ),
            "vega": Package(
                url=os.path.join(sys.prefix, "share/vega/vega.js"),
                license_url=os.path.join(sys.prefix, "share/vega/LICENSE"),
            ),
            "vega-lite": Package(
                url=os.path.join(sys.prefix, "share/vega-lite/vega-lite.js"),
                license_url=os.path.join(sys.prefix, "share/vega-lite/LICENSE"),
            ),
            "vega-embed": Package(
                url=os.path.join(sys.prefix, "share/vega-embed/vega-embed.js"),
                license_url=os.path.join(sys.prefix, "share/vega-embed/LICENSE"),
            ),
            "heroicons": Package(
                license_url=os.path.join(sys.prefix, "share/heroicons/LICENSE"),
            ),
            "prop-types": Package(
                url=os.path.join(sys.prefix, "share/prop-types/prop-types.min.js"),
                license_url=os.path.join(sys.prefix, "share/prop-types/LICENSE"),
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
