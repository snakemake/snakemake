from snakemake.ioutils.branch import branch
from snakemake.ioutils.collect import collect
from snakemake.ioutils.evaluate import evaluate
from snakemake.ioutils.exists import exists
from snakemake.ioutils.lookup import lookup
from snakemake.ioutils.prepend_param import prepend_param
from snakemake.ioutils.rule_items_proxy import rule_item_factory
from snakemake.ioutils.subpath import subpath
from snakemake.ioutils.input import parse_input, extract_checksum, flatten
from snakemake.ioutils.choose_f import choose_file, choose_folder, choose_tmp
from snakemake.ioutils.as_py_module import as_py_module


def register_in_globals(_globals):
    _globals.update(
        {
            "lookup": lookup,
            "evaluate": evaluate,
            "branch": branch,
            "collect": collect,
            "exists": exists,
            "input": rule_item_factory("input"),
            "output": rule_item_factory("output"),
            "params": rule_item_factory("params"),
            "prepend_param": prepend_param,
            "resources": rule_item_factory("resources"),
            "threads": rule_item_factory("threads"),
            "subpath": subpath,
            "parse_input": parse_input,
            "extract_checksum": extract_checksum,
            "flatten": flatten,
            "choose_file": choose_file,
            "choose_folder": choose_folder,
            "choose_tmp": choose_tmp,
            "as_py_module": as_py_module,
        }
    )
