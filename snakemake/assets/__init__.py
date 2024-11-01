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
    sha256: Optional[str] = None

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
            url="https://raw.githubusercontent.com/tailwindlabs/tailwindcss/refs/tags/v3.0.23/LICENSE",
            sha256="60e0b68c0f35c078eef3a5d29419d0b03ff84ec1df9c3f9d6e39a519a5ae7985",
        ),
        "tailwindcss/tailwind.css": Asset(
            url="https://cdn.tailwindcss.com/3.0.23?plugins=forms@0.4.0,typography@0.5.2",
            # The tailwindcss cdn checksum is not stable. Since this is only included
            # as CSS styles, the risk is low.
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
    def base_path(cls) -> Path:
        # this is called from within snakemake, so we can use importlib.resources
        if cls._base_path is None:
            cls._base_path = importlib.resources.files("snakemake.assets") / "data"
        return cls._base_path
