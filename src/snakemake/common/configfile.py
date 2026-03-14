import collections
import json
from pathlib import Path
from snakemake_interface_common.exceptions import WorkflowError


def _load_configfile(configpath_or_obj, filetype="Config"):
    "Tries to load a configfile first as JSON, then as YAML, into a dict."
    import yaml

    if isinstance(configpath_or_obj, str) or isinstance(configpath_or_obj, Path):
        obj = open(configpath_or_obj, encoding="utf-8")
    else:
        obj = configpath_or_obj

    try:
        with obj as f:
            try:
                return json.load(f, object_pairs_hook=collections.OrderedDict)
            except ValueError:
                f.seek(0)  # try again
            try:
                import yte

                return yte.process_yaml(f, require_use_yte=True)
            except yaml.YAMLError:
                raise WorkflowError(
                    f"{filetype} file is not valid JSON or YAML. "
                    "In case of YAML, make sure to not mix "
                    "whitespace and tab indentation."
                )
    except FileNotFoundError:
        raise WorkflowError(f"{filetype} file {configpath_or_obj} not found.")


def load_configfile(configpath):
    "Loads a JSON or YAML configfile as a dict, then checks that it's a dict."
    config = _load_configfile(configpath)
    if not isinstance(config, dict):
        raise WorkflowError(
            "Config file must be given as JSON or YAML with keys at top level."
        )
    return config
