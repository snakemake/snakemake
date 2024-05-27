__authors__ = "Johannes Köster, Sven Nahnsen"
__copyright__ = "Copyright 2022, Johannes Köster, Sven Nahnsen"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


import hashlib
import json

from snakemake.jobs import Job
from snakemake import script
from snakemake import wrapper
from snakemake.exceptions import WorkflowError
from snakemake.settings.enums import DeploymentMethod

# ATTENTION: increase version number whenever the hashing algorithm below changes!
__version__ = "0.1"


class ProvenanceHashMap:
    def __init__(self):
        self._hashes = dict()

    def get_provenance_hash(self, job: Job, cache_mode: str):
        versioned_hash = hashlib.sha256()
        # Ensure that semantic version changes in this module
        versioned_hash.update(self._get_provenance_hash(job, cache_mode).encode())
        versioned_hash.update(__version__.encode())
        return versioned_hash.hexdigest()

    def _get_provenance_hash(self, job: Job, cache_mode: str):
        """
        Recursively calculate hash for the output of the given job
        and all upstream jobs in a blockchain fashion.

        This is based on an idea of Sven Nahnsen.
        Fails if job has more than one output file. The reason is that there
        is no way to generate a per-output file hash without generating the files.
        This hash, however, shall work without having to generate the files,
        just by describing all steps down to a given job.
        """
        assert (cache_mode == "omit-software") or (cache_mode == "all")

        if job in self._hashes:
            return self._hashes[job]

        workflow = job.dag.workflow
        h = hashlib.sha256()

        # Hash shell command or script.
        if job.is_shell:
            # We cannot use the formatted shell command, because it also contains threads,
            # resources, and filenames (which shall be irrelevant for the hash).
            h.update(job.rule.shellcmd.encode())
        elif job.is_script:
            _, source, _, _, _ = script.get_source(
                job.rule.script,
                job.rule.workflow.sourcecache,
                basedir=job.rule.basedir,
                wildcards=job.wildcards,
                params=job.params,
            )
            h.update(source.encode())
        elif job.is_notebook:
            _, source, _, _, _ = script.get_source(
                job.rule.notebook,
                job.rule.workflow.sourcecache,
                basedir=job.rule.basedir,
                wildcards=job.wildcards,
                params=job.params,
            )
            h.update(source.encode())
        elif job.is_wrapper:
            _, source, _, _, _ = script.get_source(
                wrapper.get_script(
                    job.rule.wrapper,
                    sourcecache=job.rule.workflow.sourcecache,
                    prefix=workflow.workflow_settings.wrapper_prefix,
                ),
                job.rule.workflow.sourcecache,
                basedir=job.rule.basedir,
                wildcards=job.wildcards,
                params=job.params,
            )
            h.update(source.encode())

        # Hash params.
        for key, value in sorted(job.params._allitems()):
            if key is not None:
                h.update(key.encode())
            # If this raises a TypeError, we cannot calculate a reliable hash.
            try:
                h.update(json.dumps(value, sort_keys=True).encode())
            except TypeError as e:
                raise WorkflowError(
                    "Rule {} cannot be cached, because params "
                    "are not JSON serializable. "
                    "Consider converting them into a suitable format "
                    "if you are sure that caching is necessary. "
                    "Otherwise, deactivate caching for this rule "
                    "by removing it from the --cache command line argument "
                    "or removing the cache: true directive from the rule itself.".format(
                        job.rule.name
                    ),
                    e,
                )

        # Hash input files that are not generated by other jobs (sorted by hash value).
        for file_hash in sorted(
            hash_file(f)
            for f in job.input
            if not any(f in depfiles for depfiles in job.dag.dependencies[job].values())
        ):
            h.update(file_hash.encode())

        # Hash used containers or conda environments.
        if cache_mode != "omit-software":
            if (
                DeploymentMethod.CONDA in workflow.deployment_settings.deployment_method
                and job.conda_env
            ):
                if (
                    DeploymentMethod.APPTAINER
                    in workflow.deployment_settings.deployment_method
                    and job.conda_env.container_img_url
                ):
                    h.update(job.conda_env.container_img_url.encode())
                h.update(job.conda_env.content)
            elif (
                DeploymentMethod.APPTAINER
                in workflow.deployment_settings.deployment_method
                and job.container_img_url
            ):
                h.update(job.container_img_url.encode())

        # Generate hashes of dependencies, and add them in a blockchain fashion (as input to the current hash, sorted by hash value).
        for dep_hash in sorted(
            self._get_provenance_hash(dep, cache_mode)
            for dep in set(job.dag.dependencies[job].keys())
        ):
            h.update(dep_hash.encode())

        provenance_hash = h.hexdigest()

        # Store for reuse.
        self._hashes[job] = provenance_hash

        return provenance_hash


def hash_file(f):
    h = hashlib.sha256()
    with open(f, "rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            h.update(byte_block)
    return h.hexdigest()
