import subprocess
import os
from urllib.parse import urlparse
import hashlib

from snakemake import conda
from snakemake.common import lazy_property
from snakemake.logging import logger


class Image:
    def __init__(self, url, dag):
        self.url = url
        self._img_dir = dag.workflow.persistence.singularity_img_path

    @lazy_property
    def hash(self):
        md5hash = hashlib.md5()
        md5hash.update(self.url.encode())
        return md5hash.hexdigest()

    def pull(self, dryrun=False):
        if dryrun:
            logger.info("Singularity image {} will be pulled.".format(self.url))
            return
        logger.debug("Singularity image location: {}".format(self.path))
        scheme = urlparse(self.url).scheme
        if not scheme or scheme == "file":
            return
        if not os.path.exists(self.path):
            logger.info("Pulling singularity image {}.".format(self.url))
            try:
                subprocess.check_output(["singularity", "pull",
                                         "--name", self.hash, self.url],
                                        cwd=self._img_dir,
                                        stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                raise WorkflowError("Failed to pull singularity image "
                                    "from {}:\n{}".format(self.url,
                                                          e.stdout.decode()))


    @property
    def path(self):
        return os.path.join(self._img_dir, self.hash) + ".img"


def shellcmd(img_path, cmd, args):
    return "singularity exec {} {} bash -c \"{}\"".format(
        args, img_path, cmd.replace("\"", r"\""))
