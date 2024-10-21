# This module handles the download of non python assets.
# It should not use any modules that are not part of the standard library because it will
# be called before the setup (and dependency deployment) of the snakemake package.

from dataclasses import dataclass
import hashlib
import importlib.resources
from pathlib import Path
from typing import Dict, Optional
import urllib.request
import urllib.error


class AssetDownloadError(Exception):
    pass


@dataclass
class Asset:
    url: str
    sha256: str

    def get_content(self) -> bytes:
        """Get and validate asset content."""

        req = urllib.request.Request(self.url, headers={"User-Agent": "snakemake"})
        err = None
        for _ in range(6):
            try:
                resp = urllib.request.urlopen(req)
                content = resp.read()
            except urllib.error.URLError as e:
                err = AssetDownloadError(f"Failed to download asset {self.url}: {e}")
                continue
            if self.sha256 != hashlib.sha256(content).hexdigest():
                err = AssetDownloadError(
                    f"Checksum mismatch when downloading asset {self.url}"
                )
                continue
            return content
        assert err is not None
        raise err


class Assets:
    _base_path: Optional[Path] = None
    spec: Dict[str, Asset] = {
        "snakemake/LICENSE.md": Asset(
            url="https://raw.githubusercontent.com/snakemake/snakemake/main/LICENSE.md",
            sha256="84a1a82b05c80637744d3fe8257235c15380efa6cc32608adf4b21f17af5d2b8",
        ),
        "pygments/LICENSE": Asset(
            url="https://raw.githubusercontent.com/pygments/pygments/master/LICENSE",
            sha256="a9d66f1d526df02e29dce73436d34e56e8632f46c275bbdffc70569e882f9f17",
        ),
        "tailwindcss/LICENSE": Asset(
            url="https://raw.githubusercontent.com/tailwindlabs/tailwindcss/refs/tags/v3.0.23/LICENSE",
            sha256="60e0b68c0f35c078eef3a5d29419d0b03ff84ec1df9c3f9d6e39a519a5ae7985",
        ),
        "tailwindcss/tailwind.css": Asset(
            url="https://cdn.tailwindcss.com/3.0.23?plugins=forms@0.4.0,typography@0.5.2",
            sha256="8a597dc918fb62e05db23a5f810327a045a62c57cfda16646075138a6ac696fa",
        ),
        "react/LICENSE": Asset(
            url="https://raw.githubusercontent.com/facebook/react/refs/tags/v18.2.0/LICENSE",
            sha256="52412d7bc7ce4157ea628bbaacb8829e0a9cb3c58f57f99176126bc8cf2bfc85",
        ),
        "react/react.production.min.js": Asset(
            url="https://cdnjs.cloudflare.com/ajax/libs/react/18.2.0/umd/react.production.min.js",
            sha256="4b4969fa4ef3594324da2c6d78ce8766fbbc2fd121fff395aedf997db0a99a06",
        ),
        "react/react-dom.production.min.js": Asset(
            url="https://cdnjs.cloudflare.com/ajax/libs/react-dom/18.2.0/umd/react-dom.production.min.js",
            sha256="21758ed084cd0e37e735722ee4f3957ea960628a29dfa6c3ce1a1d47a2d6e4f7",
        ),
        "vega/vega.js": Asset(
            url="https://cdnjs.cloudflare.com/ajax/libs/vega/5.21.0/vega.js",
            sha256="b34c43055ef5d39a093e937522955dc359fbaec6c5b0259ae2de4c9da698e9fe",
        ),
        "vega/LICENSE": Asset(
            url="https://raw.githubusercontent.com/vega/vega/refs/tags/v5.21.0/LICENSE",
            sha256="b75f7ed0af20dedadf92c52bc236161bcf0d294ff2e6e34ca76403203349f71d",
        ),
        # Begin dependencies for vega, included in vega/vega.js
        # Versions from https://github.com/vega/vega/blob/v5.21.0/yarn.lock.
		# Via vega-crossfilter, vega-encode, vega-format,
		# vega-functions, vega-geo, vega-regression, vega-scale, vega-statistics,
        # vega-time, vega-transforms, vega-view; and via d3-geo (via vega-functions,
        # vega-geo, vega-projection); and via d3-geo-projection (via vega-projection)
		#
        # This has its own dependency, internmap (1.0.1); however, while we can
        # find code from vega-crossfilter and d3-array in vega.js, there does
        # not appear to be any code from internmap after tree-shaking.
        "d3-array/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-array/refs/tags/v2.12.1/LICENSE",
            sha256="0747bebeb06b10954913149be9b9a8bdf6fad3e6fdcbd9f9524e7a94c13d2cea",
        ),
        # Via vega-format
        "d3-format/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-format/refs/tags/v2.0.0/LICENSE",
            sha256="7eea8533ea92bd8c32a901e89ecb0305890905bb12711449565ddff96891146d",
        ),
        # Via vega-format
        "d3-time-format/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-time-format/refs/tags/v3.0.0/LICENSE",
            sha256="7a3cb0e5055874e67db9aa2d5fe26de23204fa994ffbad198901ffe9c812a717",
        ),
        # Via vega-time, vega-view; and via d3-time-format (via vega-format)
        "d3-time/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-time/refs/tags/v2.1.1/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
        ),
        # Via vega-encode, vega-scale
        "d3-interpolate/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-interpolate/refs/tags/v2.0.1/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
        ),
        # Via vega-functions, vega-geo; and via d3-interpolate (via vega-encode)
        "d3-color/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-color/refs/tags/v2.0.0/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
        ),
        # Via vega-scale
        "d3-scale/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-scale/refs/tags/v3.3.0/LICENSE",
            sha256="7eea8533ea92bd8c32a901e89ecb0305890905bb12711449565ddff96891146d",
        ),
        # Via vega-expression
        # The license file included in the NPM package does not exist directly
        # in https://github.com/DefinitelyTyped/DefinitelyTyped, so we use an
        # unpkg URL to reference the contents of the NPM package instead.
        "@types-estree/LICENSE": Asset(
            url="https://unpkg.com/@types/estree@0.0.50/LICENSE",
            sha256="c2cfccb812fe482101a8f04597dfc5a9991a6b2748266c47ac91b6a5aae15383",
        ),
        # Via vega-force
        "d3-force/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-force/refs/tags/v2.1.1/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
        ),
        # Via d3-force (via vega-force)
        "d3-dispatch/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-dispatch/refs/tags/v2.0.0/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
        ),
        # Via d3-force (via vega-force)
        "d3-quadtree/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-quadtree/refs/tags/v2.0.0/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
        ),
        # Via vega-view; and via d3-force (via vega-force)
        "d3-timer/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-timer/refs/tags/v2.0.0/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
        ),
        # Via vega-functions, vega-geo, vega-projection (directly and via
        # d3-geo-projection)
        "d3-geo/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-geo/refs/tags/v2.0.2/LICENSE",
            sha256="c3c9d41a75b64963748657932167ec4c56404b8fd557d4f89c9eeda6e2fdf39a",
        ),
        # Via vega-hierarchy
        "d3-hierarchy/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-hierarchy/refs/tags/v2.0.0/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
        ),
        # Via vega-loader
        # This has its own dependencies (commander 2.20.3, iconv-lite 0.4.24,
        # rw 1.3.3), but they are all used only in the command-line tools
        # rather than in the library that is bundled in vega.js.
        "d3-dsv/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-dsv/refs/tags/v2.0.0/LICENSE",
            sha256="8b5304265ccedbd17567aa14a0dc6b9bf4013fee44964c46aa54fdf8527d6a68",
        ),
        # Via vega-loader; also from vega-geo/src/util/contours.js, "Based on
        # https://github.com/topojson/topojson-client/blob/v3.0.0/src/stitch.js"
        # This has its own dependency (commander 2.20.3), but it is used only
        # in the command-line tools rather than in the library that is bundled
        # in vega.js. Dependency node-fetch 2.6.4, via vega-loader, is not
        # used for browser builds.
        "topojson-client/LICENSE": Asset(
            url="https://raw.githubusercontent.com/topojson/topojson-client/refs/tags/v3.1.0/LICENSE",
            sha256="4c4d15b635e04e691825a76db7d33f7f2033b55669a7430011694f31e6c65999",
        ),
        # Via vega-projection
        # This has its own dependencies (commander 2.20.3, resolve 1.20.0), but
        # they are all used only in the command-line tools rather than in the
        # library that is bundled in vega.js.
        "d3-geo-projection/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-geo-projection/refs/tags/v3.0.0/LICENSE",
            sha256="4108a126a74cc35d4d5ae39ca5d0cc926fa7c8ec40e459a0c6d3481c69decbd8",
        ),
        # Via vega-scenegraph (directly and via d3-shape)
        "d3-path/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-path/refs/tags/v2.0.0/LICENSE",
            sha256="b8265327ab678f554800e71188b1cc6b1ff57522e292d2b1c0be66f27cf328b6",
        ),
        # Via vega-scenegraph
        "d3-shape/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-shape/refs/tags/v2.1.0/LICENSE",
            sha256="7eea8533ea92bd8c32a901e89ecb0305890905bb12711449565ddff96891146d",
        ),
        # Via vega-voronoi
        "d3-delaunay/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-delaunay/refs/tags/v5.2.0/LICENSE",
            sha256="d19bf35c4080eea32e098a9c8d62540ffd425611ac8dfc856a233f0357d8b02e",
        ),
        # Via d3-delaunay (via vega-voronoi)
        "delaunator/LICENSE": Asset(
            url="https://raw.githubusercontent.com/mapbox/delaunator/refs/tags/v4.0.1/LICENSE",
            sha256="27043d1a6a0e1985fde12660decbbd3b23c67de900b00609c90d4f0aa492f425",
        ),
        # End dependencies for vega, included in vega/vega.js
        # Begin copied/derived/adapted code in vega, included in vega/vega.js
        # From vega-util/src/lru-cache.js, "adapted from
        # https://github.com/dominictarr/hashlru/ (MIT License)". Version 1.0.4
        # was the current version when the function was added to vega, in
        # commit f153863f3575bff1b42294c9cb065d80afb757ff on 2020-02-28.
        "hashlru/LICENSE": Asset(
            url="https://raw.githubusercontent.com/dominictarr/hashlru/refs/tags/v1.0.4/LICENSE",
            sha256="08e4bd7a280eaaf1fbcaad9dad1fb94412477fcdd1cf81660988303297e5e1d1",
        ),
        # From
        # vega-statistics/src/regression/{poly,r-squared,pow,log,linear}.js,
        # "Adapted from d3-regression by Harry Stevens". Version 1.2.1 was the
        # current version when the functions were added to vega, in commit
        # 71610e4456a3a4145435d83f8458748ba137a2a3 2019-05-10. That release was
        # not tagged in git, so we use the corresponding commit hash.
        "d3-regression/LICENSE": Asset(
            url="https://raw.githubusercontent.com/HarryStevens/d3-regression/e23d40a663dffba14b92cb42d9989de3a32894b5/LICENSE",
            sha256="d210203f922101502894baf700b9a392e323a26e4b603ab166c57e09a6e773b5",
        ),
        # From vega-statistics/src/regression/poly.js, "Adapted from
        # d3-regression by Harry Stevens [...] which was adapted from
        # regression-js by Tom Alexander" Version 2.0.1 was the current version
        # when the functions were added to vega, in commit
        # 71610e4456a3a4145435d83f8458748ba137a2a3 2019-05-10.
        "regression/LICENSE": Asset(
            url="https://raw.githubusercontent.com/Tom-Alexander/regression-js/refs/tags/2.0.1/LICENSE",
            sha256="2f932f5cfb9b042cc6c0089ee8004b33e3746ffeab879341dbd453c150524307",
        ),
        # From vega-statistics/src/regression/loess.js, "Adapted from
        # science.js by Jason Davies". Version 1.9.3 was the current version
        # when the functions were added to vega, in commit
        # 71610e4456a3a4145435d83f8458748ba137a2a3 on 2019-05-10.
        "science/LICENSE": Asset(
            url="https://raw.githubusercontent.com/jasondavies/science.js/refs/tags/v1.9.3/LICENSE",
            sha256="3bd1fdf686ffcad175daddcb36ee28599ac8f090b6cec2c7654118c8a6f3d4c9",
        ),
        # From src/quickselect.js in d3-array, "Based on
        # https://github.com/mourner/quickselect". Version 2.0.0 was the
        # current version when the file was added to d3-array, in commit
        # d447c2a31cd6aacf54a40b22c29620f7e17bbd7e on 2018-11-10.
        "quickselect/LICENSE": Asset(
            url="https://raw.githubusercontent.com/mourner/quickselect/refs/tags/v2.0.0/LICENSE",
            sha256="597034cb7c11c916ad407344ea99a0b08e3c443a6b4421460f1d23c411c69707",
        ),
        # From vega-statistics/src/normal.js, "Ported from Apache Commons
        # Math". Version 3.6.1 was the current version when this implementation
        # of erfinv was added to vega, in commit
        # c26050f21b0c95620f8f1a3094056716ec6c5aaa on 2019-09-26.
        "commons-math/LICENSE.txt": Asset(
            url="https://raw.githubusercontent.com/apache/commons-math/refs/tags/MATH_3_6_1/LICENSE.txt",
            sha256="64f23963615950bad9ddd31569a7f09afbbe11b81d94873ffd9b1cac6081a11d",
        ),
        "commons-math/NOTICE.txt": Asset(
            url="https://raw.githubusercontent.com/apache/commons-math/refs/tags/MATH_3_6_1/NOTICE.txt",
            sha256="5495442a32bfc2b93b4a8f2c34c5c218d16cca434aa5684fb953d9419120e3fa",
        ),
        # From vega-expression/src/parser.js, "The following expression parser
        # is based on Esprima (http://esprima.org/)." Unusually, the original
        # license text is reproduced in the source file, so we could get by
        # without a separate license file (and in fact the license file lacks
        # the necessary copyright statement present in the code) but we add it
        # for consistency. Version 2.2.0 was the last version of Esprima that
        # had a copyright statement in esprima.js attributing specific
        # individual contributors, as reproduced in the vega source.
        "esprima/LICENSE.BSD": Asset(
            url="https://raw.githubusercontent.com/jquery/esprima/refs/tags/2.2.0/LICENSE.BSD",
            sha256="0e74697a68cebdcd61502c30fe80ab7f9e341d995dcd452023654d57133534b1",
        ),
        # From vega-scenegraph/src/path/arc.js, "Copied from Inkscape svgtopdf,
        # thanks!" We cannot find the original implementation in Inkscape or
        # anywhere else, so we do not add a license file for this snippet.
        # From vega-scenegraph/src/path/parse.js, "Path parsing and rendering
        # code adapted from fabric.js -- Thanks!" This was added in commit
        # 82932143de7cef4187a34026689df12abaa25959 on 2018-12-20, and the
        # current release of fabric.js at that time was 2.4.5.
        "fabric/LICENSE": Asset(
            url="https://raw.githubusercontent.com/fabricjs/fabric.js/refs/tags/v2.4.5/LICENSE",
            sha256="9f6c2cc99aa9c618df93fed7d1cf7279d4e329d92dd2ce5e96173c73ce305055",
        ),
        # From vega-geo/src/util/{contours,density2D}.js, "Implementation
        # adapted from d3/d3-contour. Thanks!" When these routines were added
        # in 0ab6b730a7e576d33d00e12063855bb132194191 on 2019-11-11, the latest
        # version was 1.3.2.
        "d3-contour/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-contour/refs/tags/v1.3.2/LICENSE",
            sha256="5f5dcce265668080a60fbdc513f6f8ef21466780bcaa331e64ee39df19e63b30",
        ),
        # From vega-geo/src/util/contours.js, "Based on
        # https://github.com/mbostock/shapefile/blob/v0.6.2/shp/polygon.js"
        "shapefile/LICENSE": Asset(
            url="https://raw.githubusercontent.com/mbostock/shapefile/refs/tags/v0.6.2/LICENSE.txt",
            sha256="c16529a9d5b8802982abd714a6823344e24b0cb5131596bc343927ead605d708",
        ),
        # End copied/derived/adapted code in vega, included in vega/vega.js
        "vega-lite/vega-lite.js": Asset(
            url="https://cdnjs.cloudflare.com/ajax/libs/vega-lite/5.2.0/vega-lite.js",
            sha256="6eb7f93121cd9f44cf8640244f87c5e143f87c7a0b6cd113da4a9e41e3adf0aa",
        ),
        "vega-lite/LICENSE": Asset(
            url="https://raw.githubusercontent.com/vega/vega-lite/refs/tags/v5.2.0/LICENSE",
            sha256="f618900fd0d64046963b29f40590cdd1e341a2f41449f99110d82fd81fea808c",
        ),
        "vega-embed/vega-embed.js": Asset(
            url="https://cdnjs.cloudflare.com/ajax/libs/vega-embed/6.20.8/vega-embed.js",
            sha256="4e546c1f86eb200333606440e92f76e2940b905757018d9672cd1708e4e6ff0a",
        ),
        "vega-embed/LICENSE": Asset(
            url="https://raw.githubusercontent.com/vega/vega-embed/refs/tags/v6.20.8/LICENSE",
            sha256="32df67148f0fc3db0eb9e263a7b75d07f1eb14c61955005a4a39c6918d10d137",
        ),
        "heroicons/LICENSE": Asset(
            url="https://raw.githubusercontent.com/tailwindlabs/heroicons/refs/tags/v1.0.3/LICENSE",
            sha256="75523ddd65d9620bea09f84e89d0c373b4205a3708b8a1e9f9598a5438a3e641",
        ),
        "prop-types/prop-types.min.js": Asset(
            url="https://cdnjs.cloudflare.com/ajax/libs/prop-types/15.7.2/prop-types.min.js",
            sha256="4c88350517ee82aa4f3368e67ef1a453ca6636dcfa6449b4e3d6faa5c877066e",
        ),
        "prop-types/LICENSE": Asset(
            url="https://raw.githubusercontent.com/facebook/prop-types/refs/tags/v15.7.2/LICENSE",
            sha256="f657f99d3fb9647db92628e96007aabb46e5f04f33e49999075aab8e250ca7ce",
        ),
    }

    @classmethod
    def deploy(cls) -> None:
        # this has to work from setup.py without being able to load the snakemake
        # modules.
        base_path = Path(__file__).parent / "data"
        for asset_path, asset in cls.spec.items():
            target_path = base_path / asset_path

            if target_path.exists():
                with open(target_path, "rb") as fin:
                    # file is already present, check if it is up to date
                    if asset.sha256 == hashlib.sha256(fin.read()).hexdigest():
                        continue

            target_path.parent.mkdir(parents=True, exist_ok=True)
            with open(target_path, "wb") as fout:
                fout.write(asset.get_content())

    @classmethod
    def get_content(cls, asset_path: str) -> str:
        try:
            return (cls.base_path() / asset_path).read_text(encoding="utf-8")
        except FileNotFoundError:
            from snakemake.logging import logger

            logger.warning(
                f"Asset {asset_path} not found (development setup?), downloading..."
            )
            return cls.spec[asset_path].get_content().decode("utf-8")

    @classmethod
    def base_path(cls) -> Path:
        # this is called from within snakemake, so we can use importlib.resources
        if cls._base_path is None:
            cls._base_path = importlib.resources.files("snakemake.assets") / "data"
        return cls._base_path
