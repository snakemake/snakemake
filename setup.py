from pathlib import Path
import sys
import tomli

from setuptools import setup

source_dir = Path(__file__).parent

# import Assets, while avoiding that the rest of snakemake is imported here before
# setup has been called.
sys.path.append(str(source_dir / "snakemake"))
from assets import Assets  # type: ignore # noqa: E402

# download online assets
Assets.deploy()

# Read version from pyproject.toml
with open(source_dir / "pyproject.toml", "rb") as f:
    pyproject = tomli.load(f)
    version = pyproject.get("project").get("version")

setup(
    name="snakemake",
    version=version,
    package_data={
        "snakemake": [
            "assets/data/**/*",
            "report/html_reporter/template/**/*",
            "report/html_reporter/template/*",
            "common/tests/testcases/**/*",
        ]
    },
)
