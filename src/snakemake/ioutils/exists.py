import inspect
import os

from snakemake.common import async_run
import snakemake.io
from snakemake_interface_common.exceptions import WorkflowError


def exists(path):
    """Return True if the given file or directory exists.

    This function considers any storage arguments given to Snakemake.
    """
    func_context = inspect.currentframe().f_back.f_locals
    func_context_global = inspect.currentframe().f_back.f_globals

    workflow = func_context.get("workflow") or func_context_global.get("workflow")

    if workflow is None:
        raise WorkflowError(
            "The exists function can only be used within a Snakemake workflow "
            "(the global variable 'workflow' has to be present)."
        )

    path = workflow.modifier.path_modifier.apply_default_storage(path)
    if snakemake.io.is_flagged(path, "storage_object"):
        return async_run(path.flags["storage_object"].managed_exists())
    else:
        return os.path.exists(path)
