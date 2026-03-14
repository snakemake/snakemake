from pathlib import Path
import sys

from setuptools import setup

# ensure the current directory is on sys.path so versioneer can be imported
# when pip uses PEP 517/518 build rules.
# https://github.com/python-versioneer/python-versioneer/issues/193

source_dir = Path(__file__).parent
sys.path.append(str(source_dir))

# import Assets, while avoiding that the rest of snakemake is imported here before
# setup has been called.
sys.path.append(str(source_dir / "src" / "snakemake"))
from assets import Assets

# download online assets
Assets.deploy()

setup(
    name="snakemake",
    use_scm_version=True,
    package_data={
        "snakemake": [
            "assets/data/**/*",
            "report/html_reporter/template/**/*",
            "report/html_reporter/template/*",
            "common/tests/testcases/**/*",
            "unit_tests/templates/*",
        ]
    },
)
