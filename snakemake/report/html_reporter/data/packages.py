import json
from pathlib import Path
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

    # packages declared here must be downloaded to the prefix share/snakemake/assets
    # via setuptools_download in setup.cfg
    return Packages(
        {
            "snakemake": Package(
                version=snakemake.__version__.split("+")[0],
                license_path="snakemake/LICENSE.md",
            ),
            "pygments": Package(
                version=pygments.__version__,
                license_path="pygments/LICENSE",
            ),
            "tailwindcss": Package(
                license_path="tailwindcss/LICENSE",
                source_path="tailwindcss/tailwind.css",
            ),
            "react": Package(
                license_path="react/LICENSE",
                main="react/react.production.min.js",
                dom="react/react-dom.production.min.js",
            ),
            "vega": Package(
                source_path="vega/vega.js",
                license_path="vega/LICENSE",
            ),
            "vega-lite": Package(
                source_path="vega-lite/vega-lite.js",
                license_path="vega-lite/LICENSE",
            ),
            "vega-embed": Package(
                source_path="vega-embed/vega-embed.js",
                license_path="vega-embed/LICENSE",
            ),
            "heroicons": Package(
                license_path="heroicons/LICENSE",
            ),
            "prop-types": Package(
                source_path="prop-types/prop-types.min.js",
                license_path="prop-types/LICENSE",
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
    def __init__(
        self, version=None, license_path=None, source_path=None, **source_paths
    ):
        def read_source(path: Path):
            with path.open() as file:
                return file.read()

        self.version = version
        asset_prefix = Path(sys.prefix) / "share" / "snakemake" / "assets"
        
        self.license = read_source(asset_prefix / license_path)
        if source_path is not None:
            self.source = read_source(asset_prefix / source_path)
        else:
            self.sources = {
                name: read_source(asset_prefix / path)
                for name, path in source_paths.items()
            }

    def get_record(self):
        return {"version": self.version, "license": self.license}
