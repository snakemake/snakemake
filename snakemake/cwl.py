__author__ = "Johannes Köster"
__copyright__ = "Copyright 2018-2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from urllib.request import pathname2url
import os
import subprocess
import tempfile
import json
import shutil
import uuid
from itertools import chain

from snakemake.utils import format
from snakemake.logging import logger
from snakemake.exceptions import WorkflowError
from snakemake.shell import shell
from snakemake.common import get_container_image, Mode
from snakemake.io import is_flagged


def cwl(
    path,
    basedir,
    input,
    output,
    params,
    wildcards,
    threads,
    resources,
    log,
    config,
    rulename,
    use_singularity,
    bench_record,
    jobid,
):
    """
    Load cwl from the given basedir + path and execute it.
    """
    if shutil.which("cwltool") is None:
        raise WorkflowError(
            "'cwltool' must be in PATH in order to execute cwl directive."
        )

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


def job_to_cwl(job, dag, outputs, inputs):
    """Convert a job with its dependencies to a CWL workflow step.
    """

    if job.dynamic_output:
        raise WorkflowError("Dynamic output is not supported by CWL conversion.")
    for f in job.output:
        if os.path.isabs(f):
            raise WorkflowError(
                "All output files have to be relative to the " "working directory."
            )

    get_output_id = lambda job, i: "#main/job-{}/{}".format(job.jobid, i)

    dep_ids = {
        o: get_output_id(dep, i)
        for dep, files in dag.dependencies[job].items()
        for i, o in enumerate(dep.output)
        if o in files
    }
    files = [f for f in job.input if f not in dep_ids]
    if job.conda_env_file:
        files.add(os.path.relpath(job.conda_env_file))

    out = [get_output_id(job, i) for i, _ in enumerate(job.output)]

    def workdir_entry(i, f):
        location = "??inputs.input_files[{}].location??".format(i)
        if f.is_directory:
            entry = {
                "class": "Directory",
                "basename": os.path.basename(f),
                "location": location,
            }
        else:
            entry = {
                "class": "File",
                "basename": os.path.basename(f),
                "location": location,
            }
        return "$({})".format(
            json.dumps(outer_entry(f, entry)).replace('"??', "").replace('??"', "")
        ).replace('"', "'")

    def outer_entry(f, entry):
        parent = os.path.dirname(f)
        if parent:
            return outer_entry(
                parent,
                {
                    "class": "Directory",
                    "basename": os.path.basename(parent),
                    "listing": [entry],
                },
            )
        else:
            return entry

    if job in dag.targetjobs:
        # TODO this maps output files into the cwd after the workflow is complete.
        # We need to find a way to define subdirectories though. Otherwise,
        # there can be name clashes, and it will also become very crowded.
        outputs.append(
            {
                "type": {"type": "array", "items": "File"},
                "outputSource": "#main/job-{}/output_files".format(job.jobid),
                "id": "#main/output/job-{}".format(job.jobid),
            }
        )

    cwl = {
        "run": "#snakemake-job",
        "requirements": {
            "InitialWorkDirRequirement": {
                "listing": [
                    {"writable": True, "entry": workdir_entry(i, f)}
                    for i, f in enumerate(
                        chain(
                            files,
                            (f for dep in dag.dependencies[job] for f in dep.output),
                        )
                    )
                ]
            }
        },
        "in": {
            "cores": {"default": job.threads},
            "target_files": {"default": job.output._plainstrings()},
            "rules": {"default": [job.rule.name]},
        },
        "out": ["output_files"],
        "id": "#main/job-{}".format(job.jobid),
    }
    if files:
        inputs.append(
            {
                "type": {"type": "array", "items": "File"},
                "default": [{"class": "File", "location": f} for f in files],
                "id": "#main/input/job-{}".format(job.jobid),
            }
        )

    input_files = []
    if files:
        input_files.append("#main/input/job-{}".format(job.jobid))
    input_files.extend(
        "#main/job-{}/output_files".format(dep.jobid) for dep in dag.dependencies[job]
    )

    cwl["in"]["input_files"] = {"source": input_files, "linkMerge": "merge_flattened"}

    return cwl


def dag_to_cwl(dag):
    """Convert a given DAG to a CWL workflow, which is returned as JSON object.
    """
    snakemake_cwl = {
        "class": "CommandLineTool",
        "id": "#snakemake-job",
        "label": "Snakemake job executor",
        "hints": [{"dockerPull": get_container_image(), "class": "DockerRequirement"}],
        "baseCommand": "snakemake",
        "requirements": {"ResourceRequirement": {"coresMin": "$(inputs.cores)"}},
        "arguments": [
            "--force",
            "--keep-target-files",
            "--keep-remote",
            "--force-use-threads",
            "--wrapper-prefix",
            dag.workflow.wrapper_prefix,
            "--notemp",
            "--quiet",
            "--use-conda",
            "--no-hooks",
            "--nolock",
            "--mode",
            str(Mode.subprocess),
        ],
        "inputs": {
            "snakefile": {
                "type": "File",
                "default": {
                    "class": "File",
                    "location": os.path.relpath(dag.workflow.snakefile),
                },
                "inputBinding": {"prefix": "--snakefile"},
            },
            "sources": {
                "type": "File[]",
                "default": [
                    {"class": "File", "location": f} for f in dag.workflow.get_sources()
                ],
            },
            "cores": {
                "type": "int",
                "default": 1,
                "inputBinding": {"prefix": "--cores"},
            },
            "rules": {
                "type": "string[]?",
                "inputBinding": {"prefix": "--allowed-rules"},
            },
            "input_files": {"type": "File[]", "default": []},
            "target_files": {"type": "string[]?", "inputBinding": {"position": 0}},
        },
        "outputs": {
            "output_files": {
                "type": {"type": "array", "items": "File"},
                "outputBinding": {"glob": "$(inputs.target_files)"},
            }
        },
    }
    groups = dag.get_jobs_or_groups()
    outputs = []
    inputs = []

    dag_cwl = [job_to_cwl(job, dag, outputs, inputs) for job in groups]

    return {
        "cwlVersion": "v1.0",
        "$graph": [
            snakemake_cwl,
            {
                "class": "Workflow",
                "requirements": {
                    "InlineJavascriptRequirement": {},
                    "MultipleInputFeatureRequirement": {},
                },
                "steps": dag_cwl,
                "inputs": inputs,
                "outputs": outputs,
                "id": "#main",
            },
        ],
    }
