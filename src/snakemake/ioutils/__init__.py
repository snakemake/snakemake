from snakemake.ioutils.branch import branch
from snakemake.ioutils.collect import collect
from snakemake.ioutils.evaluate import evaluate
from snakemake.ioutils.exists import exists
from snakemake.ioutils.lookup import lookup
from snakemake.ioutils.rule_items_proxy import rule_item_factory
from snakemake.ioutils.subpath import subpath
from snakemake.ioutils.input import parse_input, extract_checksum, flatten


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
            "resources": rule_item_factory("resources"),
            "threads": rule_item_factory("threads"),
            "subpath": subpath,
            "parse_input": parse_input,
            "extract_checksum": extract_checksum,
            "flatten": flatten,
        }
    )
