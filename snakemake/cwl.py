__author__ = "Johannes Köster"
__copyright__ = "Copyright 2018, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from urllib.request import pathname2url
import os

from snakemake.utils import format
from snakemake.logging import logger
from snakemake.exceptions import WorkflowError



def cwl(path, basedir, input, output, params, wildcards, threads, resources,
        log, config, rulename, conda_env, singularity_img, singularity_args,
        bench_record):
    """
    Load a script from the given basedir + path and execute it.
    Supports Python 3 and R.
    """
    import cwltool.factory

    if not path.startswith("http"):
        if path.startswith("file://"):
            path = path[7:]
        elif path.startswith("file:"):
            path = path[5:]
        if not os.path.isabs(path):
            path = os.path.abspath(os.path.join(basedir, path))
        path = "file://" + path
    path = format(path, stepout=1)
    if path.startswith("file://"):
        sourceurl = "file:" + pathname2url(path[7:])
    else:
        sourceurl = path

    fac = cwltool.factory.Factory()

    tool = fac.make(sourceurl)

    def file_spec(f):
        return {
            "path": f,
            "class": "File",
            "contents": None,
            "basename": os.path.basename(f)
        }
    _input = dict()
    for name, f in input.items():
        _input[name] = [file_spec(f)]
    for name, p in params.items():
        _input[name] = [str(p)]

    tool(output=output[0], **_input)
