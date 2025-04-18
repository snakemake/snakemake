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
    version: Optional[str] = None
    sha256: Optional[str] = None

    def get_content(self) -> bytes:
        """Get and validate asset content."""
        url = self.url.format(version=self.version) if self.version else self.url
        req = urllib.request.Request(url, headers={"User-Agent": "snakemake"})
        err = None
        for _ in range(6):
            try:
                resp = urllib.request.urlopen(req)
                content = resp.read()
            except urllib.error.URLError as e:
                err = AssetDownloadError(f"Failed to download asset {url}: {e}")
                continue
            if self.sha256 is not None:
                content_sha = hashlib.sha256(content).hexdigest()
                if self.sha256 != content_sha:
                    err = AssetDownloadError(
                        f"Checksum mismatch when downloading asset {self.url} "
                        f"(sha: {content_sha}). First 100 bytes:\n{content[:100].decode()}"
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
            url="https://raw.githubusercontent.com/tailwindlabs/tailwindcss/refs/tags/v{version}/LICENSE",
            sha256="60e0b68c0f35c078eef3a5d29419d0b03ff84ec1df9c3f9d6e39a519a5ae7985",
            version="3.4.16",
        ),
        "tailwindcss/tailwind.css": Asset(
            url="https://cdn.tailwindcss.com/{version}?plugins=forms@0.5.9,typography@0.5.2",
            # The tailwindcss cdn checksum is not stable. Since this is only included
            # as CSS styles, the risk is low.
            version="3.4.16",
        ),
        "react/LICENSE": Asset(
            url="https://raw.githubusercontent.com/facebook/react/refs/tags/v{version}/LICENSE",
            sha256="52412d7bc7ce4157ea628bbaacb8829e0a9cb3c58f57f99176126bc8cf2bfc85",
            version="18.2.0",
        ),
        "react/react.production.min.js": Asset(
            url="https://cdnjs.cloudflare.com/ajax/libs/react/{version}/umd/react.production.min.js",
            sha256="4b4969fa4ef3594324da2c6d78ce8766fbbc2fd121fff395aedf997db0a99a06",
            version="18.2.0",
        ),
        "react/react-dom.production.min.js": Asset(
            url="https://cdnjs.cloudflare.com/ajax/libs/react-dom/{version}/umd/react-dom.production.min.js",
            sha256="21758ed084cd0e37e735722ee4f3957ea960628a29dfa6c3ce1a1d47a2d6e4f7",
            version="18.2.0",
        ),
        "vega/vega.js": Asset(
            url="https://cdnjs.cloudflare.com/ajax/libs/vega/{version}/vega.js",
            sha256="b34c43055ef5d39a093e937522955dc359fbaec6c5b0259ae2de4c9da698e9fe",
            version="5.21.0",
        ),
        "vega/LICENSE": Asset(
            url="https://raw.githubusercontent.com/vega/vega/refs/tags/v{version}/LICENSE",
            sha256="b75f7ed0af20dedadf92c52bc236161bcf0d294ff2e6e34ca76403203349f71d",
            version="5.21.0",
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
            url="https://raw.githubusercontent.com/d3/d3-array/refs/tags/v{version}/LICENSE",
            sha256="0747bebeb06b10954913149be9b9a8bdf6fad3e6fdcbd9f9524e7a94c13d2cea",
            version="2.12.1",
        ),
        # Via vega-format
        "d3-format/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-format/refs/tags/v{version}/LICENSE",
            sha256="7eea8533ea92bd8c32a901e89ecb0305890905bb12711449565ddff96891146d",
            version="2.0.0",
        ),
        # Via vega-format
        "d3-time-format/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-time-format/refs/tags/v{version}/LICENSE",
            sha256="7a3cb0e5055874e67db9aa2d5fe26de23204fa994ffbad198901ffe9c812a717",
            version="3.0.0",
        ),
        # Via vega-time, vega-view; and via d3-time-format (via vega-format)
        "d3-time/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-time/refs/tags/v{version}/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
            version="2.1.1",
        ),
        # Via vega-encode, vega-scale
        "d3-interpolate/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-interpolate/refs/tags/v{version}/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
            version="2.0.1",
        ),
        # Via vega-functions, vega-geo; and via d3-interpolate (via vega-encode)
        "d3-color/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-color/refs/tags/v{version}/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
            version="2.0.0",
        ),
        # Via vega-scale
        "d3-scale/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-scale/refs/tags/v{version}/LICENSE",
            sha256="7eea8533ea92bd8c32a901e89ecb0305890905bb12711449565ddff96891146d",
            version="3.3.0",
        ),
        # Via vega-expression (vega-lite also depends on vega-expression)
        # The license file included in the NPM package does not exist directly
        # in https://github.com/DefinitelyTyped/DefinitelyTyped, so we use an
        # unpkg URL to reference the contents of the NPM package instead.
        "@types-estree/LICENSE": Asset(
            url="https://unpkg.com/@types/estree@{version}/LICENSE",
            sha256="c2cfccb812fe482101a8f04597dfc5a9991a6b2748266c47ac91b6a5aae15383",
            version="0.0.50",
        ),
        # Via vega-force
        "d3-force/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-force/refs/tags/v{version}/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
            version="2.1.1",
        ),
        # Via d3-force (via vega-force)
        "d3-dispatch/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-dispatch/refs/tags/v{version}/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
            version="2.0.0",
        ),
        # Via d3-force (via vega-force)
        "d3-quadtree/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-quadtree/refs/tags/v{version}/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
            version="2.0.0",
        ),
        # Via vega-view; and via d3-force (via vega-force)
        "d3-timer/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-timer/refs/tags/v{version}/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
            version="2.0.0",
        ),
        # Via vega-functions, vega-geo, vega-projection (directly and via
        # d3-geo-projection)
        "d3-geo/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-geo/refs/tags/v{version}/LICENSE",
            sha256="c3c9d41a75b64963748657932167ec4c56404b8fd557d4f89c9eeda6e2fdf39a",
            version="2.0.2",
        ),
        # Via vega-hierarchy
        "d3-hierarchy/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-hierarchy/refs/tags/v{version}/LICENSE",
            sha256="e1211892da0b0e0585b7aebe8f98c1274fba15bafe47fa1f4ee8a7a502c06304",
            version="2.0.0",
        ),
        # Via vega-loader
        # This has its own dependencies (commander 2.20.3, iconv-lite 0.4.24,
        # rw 1.3.3), but they are all used only in the command-line tools
        # rather than in the library that is bundled in vega.js.
        "d3-dsv/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-dsv/refs/tags/v{version}/LICENSE",
            sha256="8b5304265ccedbd17567aa14a0dc6b9bf4013fee44964c46aa54fdf8527d6a68",
            version="2.0.0",
        ),
        # Via vega-loader; also from vega-geo/src/util/contours.js, "Based on
        # https://github.com/topojson/topojson-client/blob/v3.0.0/src/stitch.js"
        # This has its own dependency (commander 2.20.3), but it is used only
        # in the command-line tools rather than in the library that is bundled
        # in vega.js. Dependency node-fetch 2.6.4, via vega-loader, is not
        # used for browser builds.
        "topojson-client/LICENSE": Asset(
            url="https://raw.githubusercontent.com/topojson/topojson-client/refs/tags/v{version}/LICENSE",
            sha256="4c4d15b635e04e691825a76db7d33f7f2033b55669a7430011694f31e6c65999",
            version="3.1.0",
        ),
        # Via vega-projection
        # This has its own dependencies (commander 2.20.3, resolve 1.20.0), but
        # they are all used only in the command-line tools rather than in the
        # library that is bundled in vega.js.
        "d3-geo-projection/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-geo-projection/refs/tags/v{version}/LICENSE",
            sha256="4108a126a74cc35d4d5ae39ca5d0cc926fa7c8ec40e459a0c6d3481c69decbd8",
            version="3.0.0",
        ),
        # Via vega-scenegraph (directly and via d3-shape)
        "d3-path/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-path/refs/tags/v{version}/LICENSE",
            sha256="b8265327ab678f554800e71188b1cc6b1ff57522e292d2b1c0be66f27cf328b6",
            version="2.0.0",
        ),
        # Via vega-scenegraph
        "d3-shape/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-shape/refs/tags/v{version}/LICENSE",
            sha256="7eea8533ea92bd8c32a901e89ecb0305890905bb12711449565ddff96891146d",
            version="2.1.0",
        ),
        # Via vega-voronoi
        "d3-delaunay/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-delaunay/refs/tags/v{version}/LICENSE",
            sha256="d19bf35c4080eea32e098a9c8d62540ffd425611ac8dfc856a233f0357d8b02e",
            version="5.2.0",
        ),
        # Via d3-delaunay (via vega-voronoi)
        "delaunator/LICENSE": Asset(
            url="https://raw.githubusercontent.com/mapbox/delaunator/refs/tags/v{version}/LICENSE",
            sha256="27043d1a6a0e1985fde12660decbbd3b23c67de900b00609c90d4f0aa492f425",
            version="4.0.1",
        ),
        # End dependencies for vega, included in vega/vega.js
        # Begin copied/derived/adapted code in vega, included in vega/vega.js
        # From vega-util/src/lruCache.js, "adapted from
        # https://github.com/dominictarr/hashlru/ (MIT License)". Version 1.0.4
        # was the current version when the function was added to vega, in
        # commit f153863f3575bff1b42294c9cb065d80afb757ff on 2020-02-28.
        "hashlru/LICENSE": Asset(
            url="https://raw.githubusercontent.com/dominictarr/hashlru/refs/tags/v{version}/LICENSE",
            sha256="08e4bd7a280eaaf1fbcaad9dad1fb94412477fcdd1cf81660988303297e5e1d1",
            version="1.0.4",
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
            version="1.2.1",
        ),
        # From vega-statistics/src/regression/poly.js, "Adapted from
        # d3-regression by Harry Stevens [...] which was adapted from
        # regression-js by Tom Alexander" Version 2.0.1 was the current version
        # when the functions were added to vega, in commit
        # 71610e4456a3a4145435d83f8458748ba137a2a3 2019-05-10.
        "regression/LICENSE": Asset(
            url="https://raw.githubusercontent.com/Tom-Alexander/regression-js/refs/tags/{version}/LICENSE",
            sha256="2f932f5cfb9b042cc6c0089ee8004b33e3746ffeab879341dbd453c150524307",
            version="2.0.1",
        ),
        # From vega-statistics/src/regression/loess.js, "Adapted from
        # science.js by Jason Davies". Version 1.9.3 was the current version
        # when the functions were added to vega, in commit
        # 71610e4456a3a4145435d83f8458748ba137a2a3 on 2019-05-10.
        "science/LICENSE": Asset(
            url="https://raw.githubusercontent.com/jasondavies/science.js/refs/tags/v{version}/LICENSE",
            sha256="3bd1fdf686ffcad175daddcb36ee28599ac8f090b6cec2c7654118c8a6f3d4c9",
            version="1.9.3",
        ),
        # From src/quickselect.js in d3-array, "Based on
        # https://github.com/mourner/quickselect". Version 2.0.0 was the
        # current version when the file was added to d3-array, in commit
        # d447c2a31cd6aacf54a40b22c29620f7e17bbd7e on 2018-11-10.
        "quickselect/LICENSE": Asset(
            url="https://raw.githubusercontent.com/mourner/quickselect/refs/tags/v{version}/LICENSE",
            sha256="597034cb7c11c916ad407344ea99a0b08e3c443a6b4421460f1d23c411c69707",
            version="2.0.0",
        ),
        # From vega-statistics/src/normal.js, "Ported from Apache Commons
        # Math". Version 3.6.1 was the current version when this implementation
        # of erfinv was added to vega, in commit
        # c26050f21b0c95620f8f1a3094056716ec6c5aaa on 2019-09-26.
        "commons-math/LICENSE.txt": Asset(
            url="https://raw.githubusercontent.com/apache/commons-math/refs/tags/MATH_3_6_1/LICENSE.txt",
            sha256="64f23963615950bad9ddd31569a7f09afbbe11b81d94873ffd9b1cac6081a11d",
            version="3.6.1",
        ),
        "commons-math/NOTICE.txt": Asset(
            url="https://raw.githubusercontent.com/apache/commons-math/refs/tags/MATH_3_6_1/NOTICE.txt",
            sha256="5495442a32bfc2b93b4a8f2c34c5c218d16cca434aa5684fb953d9419120e3fa",
            version="3.6.1",
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
            url="https://raw.githubusercontent.com/jquery/esprima/refs/tags/{version}/LICENSE.BSD",
            sha256="0e74697a68cebdcd61502c30fe80ab7f9e341d995dcd452023654d57133534b1",
            version="2.2.0",
        ),
        # From vega-scenegraph/src/path/arc.js, "Copied from Inkscape svgtopdf,
        # thanks!" We cannot find the original implementation in Inkscape or
        # anywhere else, so we do not add a license file for this snippet.
        # From vega-scenegraph/src/path/parse.js, "Path parsing and rendering
        # code adapted from fabric.js -- Thanks!" This was added in commit
        # 82932143de7cef4187a34026689df12abaa25959 on 2018-12-20, and the
        # current release of fabric.js at that time was 2.4.5.
        "fabric/LICENSE": Asset(
            url="https://raw.githubusercontent.com/fabricjs/fabric.js/refs/tags/v{version}/LICENSE",
            sha256="9f6c2cc99aa9c618df93fed7d1cf7279d4e329d92dd2ce5e96173c73ce305055",
            version="2.4.5",
        ),
        # From vega-geo/src/util/{contours,density2D}.js, "Implementation
        # adapted from d3/d3-contour. Thanks!" When these routines were added
        # in 0ab6b730a7e576d33d00e12063855bb132194191 on 2019-11-11, the latest
        # version was 1.3.2.
        "d3-contour/LICENSE": Asset(
            url="https://raw.githubusercontent.com/d3/d3-contour/refs/tags/v{version}/LICENSE",
            sha256="5f5dcce265668080a60fbdc513f6f8ef21466780bcaa331e64ee39df19e63b30",
            version="1.3.2",
        ),
        # From vega-geo/src/util/contours.js, "Based on
        # https://github.com/mbostock/shapefile/blob/v0.6.2/shp/polygon.js"
        "shapefile/LICENSE": Asset(
            url="https://raw.githubusercontent.com/mbostock/shapefile/refs/tags/v{version}/LICENSE.txt",
            sha256="c16529a9d5b8802982abd714a6823344e24b0cb5131596bc343927ead605d708",
            version="0.6.2",
        ),
        # End copied/derived/adapted code in vega, included in vega/vega.js
        "vega-lite/vega-lite.js": Asset(
            url="https://cdnjs.cloudflare.com/ajax/libs/vega-lite/{version}/vega-lite.js",
            sha256="6eb7f93121cd9f44cf8640244f87c5e143f87c7a0b6cd113da4a9e41e3adf0aa",
            version="5.2.0",
        ),
        "vega-lite/LICENSE": Asset(
            url="https://raw.githubusercontent.com/vega/vega-lite/refs/tags/v{version}/LICENSE",
            sha256="f618900fd0d64046963b29f40590cdd1e341a2f41449f99110d82fd81fea808c",
            version="5.2.0",
        ),
        # Begin dependencies for vega-lite, included in vega-lite/vega-lite.js
        # Versions from https://github.com/vega/vega-lite/blob/v5.2.0/yarn.lock.
        "@types-clone/LICENSE": Asset(
            url="https://unpkg.com/@types/clone@{version}/LICENSE",
            sha256="c2cfccb812fe482101a8f04597dfc5a9991a6b2748266c47ac91b6a5aae15383",
            version="2.1.1",
        ),
        "array-flat-polyfill/LICENSE": Asset(
            # Releases are not tagged in git; we use the commit hash
            # corresponding to the 1.0.1 release
            url="https://raw.githubusercontent.com/jonathantneal/array-flat-polyfill/362d855cb6ea2ef12f0676b116c5f9e4233b6f31/LICENSE.md",
            sha256="597756adcb51f243ef4fb386920377f61d012ace0904364e1a8ee9aaec6afc84",
            version="1.0.1",
        ),
        "clone/LICENSE": Asset(
            url="https://raw.githubusercontent.com/pvorb/clone/refs/tags/v{version}/LICENSE",
            sha256="3fb0857ef0133928cf72c88dfc464e931486e88778961eedec25585e2321507f",
            version="2.1.2",
        ),
        "fast-deep-equal/LICENSE": Asset(
            url="https://raw.githubusercontent.com/epoberezkin/fast-deep-equal/refs/tags/v{version}/LICENSE",
            sha256="7bf9b2de73a6b356761c948d0e9eeb4be6c1270bd04c79cd489c1e400ffdfc1a",
            version="3.1.3",
        ),
        "fast-json-stable-stringify/LICENSE": Asset(
            url="https://raw.githubusercontent.com/epoberezkin/fast-json-stable-stringify/refs/tags/v{version}/LICENSE",
            sha256="a833d366242c298cf1b10263516572fb8dcbe68eb5072cdcac2b4546e2b4eb36",
            version="2.1.0",
        ),
        "json-stringify-pretty-compact/LICENSE": Asset(
            url="https://raw.githubusercontent.com/lydell/json-stringify-pretty-compact/refs/tags/v{version}/LICENSE",
            sha256="d0800a7c9b5f723f7fe029f75d3a9ae44173178b9729d1521e18371ab24ae97a",
            version="3.0.0",
        ),
        # The tslib package is among vega-lite's direct dependencies, and the
        # yarn.lock file shows three versions could be included: 2.2.0, 2.1.0
        # (both MIT-licensed), and 1.1.0 (Apache-2.0-licensed). However, we
        # found no trace of any of these in the actual vega-lite.js bundle.
        # ----
        # The following packages are part of vega, and are already covered by
        # vega/LICENSE.
        # - vega-event-selector: no dependencies
        # - vega-expression: depends on vega-util, @types/estree (already
        #   covered as one of vega's dependencies)
        # - vega-util (no dependencies; includes code from hashlru, already
        #   covered as one of vega's dependencies)
        # ----
        # End dependencies for vega-lite, included in vega-lite/vega-lite.js
        # Begin copied/derived/adapted code in vega-lite, included in vega-lite/vega-lite.js
        # - hashlru is included in vega-util; already covered as one of vega's dependencies
        # End copied/derived/adapted code in vega-lite, included in vega-lite/vega-lite.js
        "vega-embed/vega-embed.js": Asset(
            url="https://cdnjs.cloudflare.com/ajax/libs/vega-embed/{version}/vega-embed.js",
            sha256="4e546c1f86eb200333606440e92f76e2940b905757018d9672cd1708e4e6ff0a",
            version="6.20.8",
        ),
        "vega-embed/LICENSE": Asset(
            url="https://raw.githubusercontent.com/vega/vega-embed/refs/tags/v{version}/LICENSE",
            sha256="32df67148f0fc3db0eb9e263a7b75d07f1eb14c61955005a4a39c6918d10d137",
            version="6.20.8",
        ),
        # Begin dependencies for vega-embed, included in vega-embed/vega-embed.js
        # Versions from https://github.com/vega/vega-embed/blob/v6.20.8/yarn.lock.
        "fast-json-patch/LICENSE": Asset(
            url="https://raw.githubusercontent.com/Starcounter-Jack/JSON-Patch/refs/tags/v{version}/LICENSE.txt",
            sha256="26593d78926902a82eb9bc7e40d1811fc381055f23516b92fe48ca8399bf02df",
            version="3.1.0",
        ),
        # json-stringify-pretty-compact 3.0.0 is already covered as one of
        # vega's dependencies
        # ----
        # There are four versions of the semver package in yarn.lock (5.7.1,
        # 6.3.0, 7.0.0, 7.3.5), but it appears that only one copy remains in
        # the vega-embed.js bundle. The bundle has a compareBuild function,
        # added in 6.1.0. It does not export a "tokens" list, added in 6.3.0
        # and present in 7.0.0, but that is not conclusive because the package
        # was split into separate source files to improve tree-shaking. By
        # 7.3.5, the first line of the SemVer class constructor changes to
        # "options = parseOptions(options)", and this does match what we see in
        # the bundle, so we conclude that we have semver 7.3.5. In any case,
        # the ISC license was unchanged across all four versions.
        "semver/LICENSE": Asset(
            url="https://raw.githubusercontent.com/npm/node-semver/refs/tags/v{version}/LICENSE",
            sha256="4ec3d4c66cd87f5c8d8ad911b10f99bf27cb00cdfcff82621956e379186b016b",
            version="7.3.5",
        ),
        # The tslib package is among vega-embed's direct dependencies, and the
        # yarn.lock file shows four versions could be included: 2.3.1, 2.1.0,
        # 1.14.1 (all MIT-licensed), and 1.1.0 (Apache-2.0-licensed). However,
        # we found no trace of any of these in the actual vega-embed.js bundle.
        # ----
        # The following packages are part of vega, and are already covered by
        # vega/LICENSE.
        # - vega-interpreter: no dependencies
        # ----
        # This is not part of https://github.com/vega/vega, as far as we can
        # tell, so we handle it separately. The repository
        # https://github.com/vega/schema linked from NPM appears to contain
        # only JSON schemas, not JavaScript source code, so we use an unpkg URL
        # to reference the contents of the NPM package instead.
        "vega-schema-url-parser/LICENSE": Asset(
            url="https://unpkg.com/vega-schema-url-parser@{version}/LICENSE",
            sha256="eeb1606fda238f623c36d1eb5e27c33b986a446445340b220c141097424af045",
            version="2.2.0",
        ),
        "vega-themes/LICENSE": Asset(
            url="https://raw.githubusercontent.com/vega/vega-themes/refs/tags/v{version}/LICENSE",
            sha256="ceb75fa4fbebbc381c0465442a502ffc13b73f36ccb385a66bcd5f55869e3979",
            version="2.10.0",
        ),
        # Depends on vega-util, which is already covered in vega
        "vega-tooltip/LICENSE": Asset(
            url="https://raw.githubusercontent.com/vega/vega-tooltip/refs/tags/v{version}/LICENSE",
            sha256="93cd20d9cee8c77c080d1a585e3b07644ad2c8e9ccdfc86e350136884d62d935",
            version="0.28.0",
        ),
        # End dependencies for vega-embed, included in vega-embed/vega-embed.js
        # Begin copied/derived/adapted code in vega-embed, included in vega-embed/vega-embed.js
        # - hashlru is included in vega-util (via vega-tooltip); already
        #   covered as one of vega's dependencies
        # - _areEquals() is based on fast-deep-equal; already covered as one of
        #   vega-lite's dependencies
        # - googlecharts theme in vega-themes carries its own copyright statement,
        #     Copyright 2020 Google LLC.
        #
        #     Use of this source code is governed by a BSD-style
        #     license that can be found in the LICENSE file or at
        #     https://developers.google.com/open-source/licenses/bsd
        #   No plain-text version of that URL is available, and the license has
        #   the same BSD-3-Clause permission statement as vega-themes itself,
        #   so we consider the combination of vega-themes/LICENSE and the
        #   copyright statement in the source code sufficient.
        # End copied/derived/adapted code in vega-embed, included in vega-embed/vega-embed.js
        "heroicons/LICENSE": Asset(
            url="https://raw.githubusercontent.com/tailwindlabs/heroicons/refs/tags/v{version}/LICENSE",
            sha256="75523ddd65d9620bea09f84e89d0c373b4205a3708b8a1e9f9598a5438a3e641",
            version="1.0.3",
        ),
        "prop-types/prop-types.min.js": Asset(
            url="https://cdnjs.cloudflare.com/ajax/libs/prop-types/{version}/prop-types.min.js",
            sha256="4c88350517ee82aa4f3368e67ef1a453ca6636dcfa6449b4e3d6faa5c877066e",
            version="15.7.2",
        ),
        "prop-types/LICENSE": Asset(
            url="https://raw.githubusercontent.com/facebook/prop-types/refs/tags/v{version}/LICENSE",
            sha256="f657f99d3fb9647db92628e96007aabb46e5f04f33e49999075aab8e250ca7ce",
            version="15.7.2",
        ),
        # Begin dependencies for prop-types, included in prop-types/prop-types.min.js
        # Versions from https://github.com/facebook/prop-types/blob/v15.7.2/yarn.lock
        # ----
        # The loose-envify package is among prop-types's direct dependencies,
        # and the yarn.lock file shows two versions could be included: 1.4.0
        # and 1.3.1 (both MIT-licensed). However, we found no trace of
        # loose-envify code in the actual prop-types.min.js bundle, so it
        # appears this dependency was removed in tree-shaking, along with its
        # dependency js-tokens (3.0.0 and/or 4.0.0).
        # ----
        # The object-assign package (version 4.1.1) is among prop-types's
        # direct dependencies, but we found no trace of it in the actual
        # prop-types.min.js bundle, so it appears this dependency (a polyfill)
        # was removed in tree-shaking.
        # ----
        # The react-is package (version 16.8.1) is among prop-types's
        # direct dependencies, but we found no trace of it in the actual
        # prop-types.min.js bundle, so it appears this dependency was removed
        # in tree-shaking.
        # ----
        # End dependencies for prop-types, included in prop-types/prop-types.min.js
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
                    if (asset.sha256 is not None) and (
                        asset.sha256 == hashlib.sha256(fin.read()).hexdigest()
                    ):
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
    def get_version(cls, asset_path: str) -> Optional[str]:
        if asset_path in cls.spec:
            return cls.spec[asset_path].version
        else:
            return None

    @classmethod
    def base_path(cls) -> Path:
        # this is called from within snakemake, so we can use importlib.resources
        if cls._base_path is None:
            cls._base_path = importlib.resources.files("snakemake.assets") / "data"
        return cls._base_path
