import json
import snakemake

from snakemake.assets import AssetDownloadError, Assets
from snakemake_interface_common.exceptions import WorkflowError


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
            # Begin dependencies for vega, included in vega/vega.js
            "d3-array": Package(
                license_path="d3-array/LICENSE",
            ),
            "d3-format": Package(
                license_path="d3-format/LICENSE",
            ),
            "d3-time-format": Package(
                license_path="d3-time-format/LICENSE",
            ),
            "d3-time": Package(
                license_path="d3-time/LICENSE",
            ),
            "d3-interpolate": Package(
                license_path="d3-interpolate/LICENSE",
            ),
            "d3-color": Package(
                license_path="d3-color/LICENSE",
            ),
            "d3-scale": Package(
                license_path="d3-scale/LICENSE",
            ),
            "@types-estree": Package(
                license_path="@types-estree/LICENSE",
            ),
            "d3-force": Package(
                license_path="d3-force/LICENSE",
            ),
            "d3-dispatch": Package(
                license_path="d3-dispatch/LICENSE",
            ),
            "d3-quadtree": Package(
                license_path="d3-quadtree/LICENSE",
            ),
            "d3-timer": Package(
                license_path="d3-timer/LICENSE",
            ),
            "d3-geo": Package(
                license_path="d3-geo/LICENSE",
            ),
            "d3-hierarchy": Package(
                license_path="d3-hierarchy/LICENSE",
            ),
            "d3-dsv": Package(
                license_path="d3-dsv/LICENSE",
            ),
            "topojson-client": Package(
                license_path="topojson-client/LICENSE",
            ),
            "d3-geo-projection": Package(
                license_path="d3-geo-projection/LICENSE",
            ),
            "d3-path": Package(
                license_path="d3-path/LICENSE",
            ),
            "d3-shape": Package(
                license_path="d3-shape/LICENSE",
            ),
            "d3-delaunay": Package(
                license_path="d3-delaunay/LICENSE",
            ),
            "delaunator": Package(
                license_path="delaunator/LICENSE",
            ),
            # End dependencies for vega, included in vega/vega.js
            # Begin copied/derived/adapted code in vega, included in vega/vega.js
            "hashlru": Package(
                license_path="hashlru/LICENSE",
            ),
            "d3-regression": Package(
                license_path="d3-regression/LICENSE",
            ),
            "regression": Package(
                license_path="regression/LICENSE",
            ),
            "science": Package(
                license_path="science/LICENSE",
            ),
            "quickselect": Package(
                license_path="quickselect/LICENSE",
            ),
            "commons-math": Package(
                license_path="commons-math/LICENSE.txt",
                # This must be included with LICENSE.txt.
                notice="commons-math/NOTICE.txt",
            ),
            "esprima": Package(
                license_path="esprima/LICENSE.BSD",
            ),
            "fabric": Package(
                license_path="fabric/LICENSE",
            ),
            "d3-contour": Package(
                license_path="d3-contour/LICENSE",
            ),
            # End copied/derived/adapted code in vega, included in vega/vega.js
            "vega-lite": Package(
                source_path="vega-lite/vega-lite.js",
                license_path="vega-lite/LICENSE",
            ),
            # Begin dependencies for vega-lite, included in vega-lite/vega-lite.js
            # (excluding those shared with vega and therefore already documented)
            "@types-clone": Package(
                license_path="@types-clone/LICENSE",
            ),
            "array-flat-polyfill": Package(
                license_path="array-flat-polyfill/LICENSE",
            ),
            "clone": Package(
                license_path="clone/LICENSE",
            ),
            "fast-deep-equal": Package(
                license_path="fast-deep-equal/LICENSE",
            ),
            "fast-json-stable-stringify": Package(
                license_path="fast-json-stable-stringify/LICENSE",
            ),
            "json-stringify-pretty-compact": Package(
                license_path="json-stringify-pretty-compact/LICENSE",
            ),
            # End dependencies for vega-lite, included in vega-lite/vega-lite.js
            # Any copied/derived/adapted code in vega-lite is shared with vega.
            "vega-embed": Package(
                source_path="vega-embed/vega-embed.js",
                license_path="vega-embed/LICENSE",
            ),
            # Begin dependencies for vega-embed, included in vega-embed/vega-embed.js
            # (excluding those shared with vega/vega-lite and therefore already documented)
            "fast-json-patch": Package(
                license_path="fast-json-patch/LICENSE",
            ),
            "semver": Package(
                license_path="semver/LICENSE",
            ),
            "vega-schema-url-parser": Package(
                license_path="vega-schema-url-parser/LICENSE",
            ),
            "vega-themes": Package(
                license_path="vega-themes/LICENSE",
            ),
            "vega-tooltip": Package(
                license_path="vega-tooltip/LICENSE",
            ),
            # End dependencies for vega-embed, included in vega-embed/vega-embed.js
            # Any copied/derived/adapted code in vega-embed is shared with vega and/or vega-lite.
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
        self.version = (
            version
            or Assets.get_version(source_path)
            or Assets.get_version(license_path)
        )

        try:
            self.license = Assets.get_content(license_path)
            if source_path is not None:
                self.source = Assets.get_content(source_path)
            else:
                self.sources = {
                    name: Assets.get_content(path)
                    for name, path in source_paths.items()
                }
        except AssetDownloadError as e:
            raise WorkflowError(e)

    def get_record(self):
        return {"version": self.version, "license": self.license}
