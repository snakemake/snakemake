import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from common import *


# run tibanna test before any moto-related tests because they apparently render AWS environment variables invalid or uneffective.
def test_tibanna():
    workdir = dpath("test_tibanna")
    subprocess.check_call(["python", "cleanup.py"], cwd=workdir)

    os.environ["TEST_ENVVAR1"] = "test"
    os.environ["TEST_ENVVAR2"] = "test"

    run(
        workdir,
        use_conda=True,
        default_remote_prefix="snakemake-tibanna-test/1",
        tibanna=True,
        tibanna_sfn="tibanna_unicorn_johannes",
        tibanna_config=["spot_instance=true"],
    )
