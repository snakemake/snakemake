__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import inspect
import os

from snakemake.utils import format


class REncoder:
    """Encoding Pyton data structures into R."""

    @classmethod
    def encode_value(cls, value):
        if isinstance(value, list):
            return cls.encode_list(value)
        elif isinstance(value, dict):
            return cls.encode_dict(value)
        elif isinstance(value, str):
            return '"{}"'.format(value)
        elif isinstance(value, int) or isinstance(value, bool) or isinstance(
            value, float):
            return str(value)
        else:
            raise ValueError(
                "Unsupported value for conversion into R: {}".format(value))

    @classmethod
    def encode_list(cls, l):
        if not l:
            return ""
        return "c({})".format(", ".join(map(cls.encode_value, l)))

    @classmethod
    def encode_items(cls, items):
        def encode_item(item):
            name, value = item
            return "{} = {}".format(name, cls.encode_value(value))

        return ", ".join(map(encode_item, items))

    @classmethod
    def encode_dict(cls, d):
        return "list({})".format(cls.encode_items(d.items()))

    @classmethod
    def encode_namedlist(cls, namedlist):
        positional = cls.encode_list(namedlist)
        named = cls.encode_items(namedlist.items())
        source = "list("
        if positional:
            source += positional
        if named:
            source += ", " + named
        source += ")"
        return source


class Snakemake:
    def __init__(self, input, output, params, wildcards, threads, resources,
                 log, config):
        self.input = input
        self.output = output
        self.params = params
        self.wildcards = wildcards
        self.threads = threads
        self.resources = resources
        self.log = log
        self.config = config


def script(basedir, path, input, output, params, wildcards, threads, resources,
           log, config):
    """
    Load a script from the given basedir + path and execute it.
    Supports Python 3 and R.
    """
    path = format(os.path.join(basedir, path), stepout=1)

    if path.endswith(".py"):
        with open(path) as source:
            exec(compile(source.read(), path, "exec"), {
                "snakemake": Snakemake(input, output, params, wildcards,
                                       threads, resources, log, config)
            })
    elif path.endswith(".R"):
        try:
            import rpy2.robjects as robjects
        except ImportError:
            raise ValueError(
                "Python 3 package rpy2 needs to be installed to use the R function.")
        with open(path) as source:
            preamble = """
            Snakemake <- setClass(
                "Snakemake",
                slots = c(
                    input = "list",
                    output = "list",
                    params = "list",
                    wildcards = "list",
                    threads = "numeric",
                    log = "list",
                    resources = "list",
                    config = "list"
                )
            )
            snakemake <- Snakemake(
                input = {},
                output = {},
                params = {},
                wildcards = {},
                threads = {},
                log = {},
                resources = {},
                config = {}
            )
            """.format(REncoder.encode_namedlist(input),
                       REncoder.encode_namedlist(output),
                       REncoder.encode_namedlist(params),
                       REncoder.encode_namedlist(wildcards), threads,
                       REncoder.encode_namedlist(log),
                       REncoder.encode_namedlist({
                           name: value
                           for name, value in resources.items()
                           if name != "_cores" and name != "_nodes"
                       }), REncoder.encode_dict(config))
            source = preamble + source.read()
            robjects.r(source)
    else:
        raise ValueError(
            "Unsupported script: Expecting either Python (.py) or R (.R) script.")
