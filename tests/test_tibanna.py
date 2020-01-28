import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from common import *

# run tibanna test before any moto-related tests because they apparently render AWS environment variables invalid or uneffective.
def test_tibanna():
    workdir = dpath("test_tibanna")
    subprocess.check_call(["python", "cleanup.py"], cwd=workdir)
    run(
        workdir,
        use_conda=True,
        configfiles=[os.path.join(workdir, "config.json")],
        default_remote_prefix="snakemake-tibanna-test/1",
        tibanna_sfn="tibanna_unicorn_johannes",
    )
