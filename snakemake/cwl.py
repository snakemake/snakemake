__author__ = "Johannes Köster"
__copyright__ = "Copyright 2018, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from urllib.request import pathname2url
import os
import subprocess
import tempfile
import json
import shutil

from snakemake.utils import format
from snakemake.logging import logger
from snakemake.exceptions import WorkflowError
from snakemake.shell import shell


def cwl(path, basedir, input, output, params, wildcards, threads, resources,
        log, config, rulename, use_singularity, bench_record, jobid):
    """
    Load cwl from the given basedir + path and execute it.
    """
    if shutil.which("cwltool") is None:
        raise WorkflowError("'cwltool' must be in PATH in order to execute cwl directive.")

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

    def file_spec(f):
        if isinstance(f, str):
            return {"path": os.path.abspath(f), "class": "File"}
        return [file_spec(f_) for f_ in f]

    inputs = dict()
    inputs.update({name: file_spec(f) for name, f in input.items()})
    inputs.update({name: p for name, p in params.items()})
    inputs.update({name: f for name, f in output.items()})
    inputs.update({name: f for name, f in log.items()})

    args = "--singularity" if use_singularity else ""

    with tempfile.NamedTemporaryFile(mode="w") as input_file:
        json.dump(inputs, input_file)
        input_file.flush()
        cmd = "cwltool {} {} {}".format(args, sourceurl, input_file.name)
        shell(cmd, bench_record=bench_record)
