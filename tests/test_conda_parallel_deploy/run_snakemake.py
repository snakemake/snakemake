"""Run Snakemake, and hold each deployment until the test releases it.

Used by the parallel deployment tests in tests_using_conda.py. A real deployment
can be too fast for the test to see it, thus this script pauses after each
deployment, while Snakemake holds the deployment lock and before it writes the
done flag:

1. Touch the file in SNAKEMAKE_TEST_DEPLOY_STARTED.
2. Wait until the file in SNAKEMAKE_TEST_DEPLOY_RELEASE exists.

If the variables are not set, the deployment does not pause.
"""

import asyncio
import os
import time
from pathlib import Path

from snakemake_interface_software_deployment_plugins import DeployableEnvBase

from snakemake.cli import main

RELEASE_TIMEOUT = 900

_managed_deploy = DeployableEnvBase.managed_deploy


async def _paused_managed_deploy(self) -> None:
    await _managed_deploy(self)
    started = os.environ.get("SNAKEMAKE_TEST_DEPLOY_STARTED")
    release = os.environ.get("SNAKEMAKE_TEST_DEPLOY_RELEASE")
    if started:
        Path(started).touch()
    if release:
        deadline = time.time() + RELEASE_TIMEOUT
        while not Path(release).exists():
            if time.time() > deadline:
                raise TimeoutError(f"The test did not release {release}.")
            await asyncio.sleep(0.1)


DeployableEnvBase.managed_deploy = _paused_managed_deploy

if __name__ == "__main__":
    main()
