from typing import Optional
from typing import Dict
from snakemake.settings.types import NotebookEditMode
from typing import TYPE_CHECKING

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import collections
import itertools
import json
import math
import os
import pickle
import re
import shlex
import sys
import tempfile
import textwrap
from abc import ABC, abstractmethod
from collections.abc import Iterable
from numbers import Complex, Integral, Number, Real
from pathlib import Path
from typing import List, Pattern, Tuple
from urllib.error import URLError

from snakemake import iocontainers
from snakemake.iocontainers import Snakemake
from snakemake.common.misc import get_snakemake_searchpaths
from snakemake.common.constants import (
    MIN_PY_VERSION,
    ON_WINDOWS,
)
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.shell import shell
from snakemake.sourcecache import (
    LocalSourceFile,
    SourceCache,
    SourceFile,
    infer_source_file,
)
from snakemake.utils import format

if TYPE_CHECKING:
    from snakemake.executors.local import RunArgs

# TODO use this to find the right place for inserting the preamble
PY_PREAMBLE_RE = re.compile(r"from( )+__future__( )+import.*?(?P<end>[;\n])")


class REncoder:
    """Encoding Python data structures into R."""

    @classmethod
    def encode_numeric(cls, value):
        if value is None:
            return "as.numeric(NA)"
        elif isinstance(value, Integral):
            return f"{value}L"
        elif isinstance(value, Real):
            if value == float("inf"):
                return "Inf"
            elif value == float("-inf"):
                return "-Inf"
            elif math.isnan(value):
                return "NaN"
            else:
                return f"{value}"
        elif isinstance(value, Complex):
            return f"{cls.encode_numeric(value.real)}+{cls.encode_numeric(value.imag)}i"
        else:
            raise ValueError("Value is not a proper number")

    @classmethod
    def encode_value(cls, value):
        if value is None:
            return "NULL"
        elif isinstance(value, bool):
            return "TRUE" if value else "FALSE"
        elif isinstance(value, (int, float, complex)):
            return cls.encode_numeric(value)
        elif isinstance(value, str):
            return repr(value)
        elif isinstance(value, Path):
            return repr(str(value))
        elif isinstance(value, dict):
            return cls.encode_dict(value)
        elif isinstance(value, collections.abc.Iterable):
            return cls.encode_list(value)
        else:
            # Try to convert from numpy if numpy is present
            try:
                import numpy as np

                if isinstance(value, np.number):
                    return cls.encode_numeric(value.item())
                elif isinstance(value, np.bool_):
                    return "TRUE" if value else "FALSE"
            except ImportError:
                pass
        raise ValueError(f"Unsupported value for conversion into R: {value}")

    @classmethod
    def encode_list(cls, l):
        """Encode as vector if the type is homogeneous, otherwise use a list."""
        is_homogeneous = False
        if len(l) == 0:
            # An empty list is always homogeneous
            is_homogeneous = True
        else:
            # Numbers of different type can be stored in the same vector,
            # casting without information loss is acceptable (e.g. int -> float)
            for vector_type in (Number, bool, str, bytes):
                if all([isinstance(e, vector_type) for e in l]):
                    is_homogeneous = True

        if is_homogeneous:
            return "c({})".format(", ".join(map(cls.encode_value, l)))
        else:
            return "list({})".format(", ".join(map(cls.encode_value, l)))

    @classmethod
    def encode_items(cls, items):
        def encode_item(item):
            name, value = item
            return f'"{name}" = {cls.encode_value(value)}'

        return ", ".join(map(encode_item, items))

    @classmethod
    def encode_dict(cls, d):
        d = f"list({cls.encode_items(d.items())})"
        return d

    @classmethod
    def encode_namedlist(cls, namedlist: iocontainers.Namedlist):
        positional = ", ".join(map(cls.encode_value, namedlist))
        named = cls.encode_items(namedlist.items())
        source = "list("
        if positional:
            source += positional
        if named:
            source += ", " + named
        source += ")"
        return source


class JuliaEncoder:
    """Encoding Python data structures into Julia."""

    @classmethod
    def encode_value(cls, value):
        if value is None:
            return "nothing"
        elif isinstance(value, str):
            # JS string quoting works OK for Julia - see tests/test_script/scripts/test.jl
            return json.dumps(value)
        elif isinstance(value, Path):
            return json.dumps(str(value))
        elif isinstance(value, dict):
            return cls.encode_dict(value)
        elif isinstance(value, bool):
            return "true" if value else "false"
        elif isinstance(value, int) or isinstance(value, float):
            return str(value)
        elif isinstance(value, collections.abc.Iterable):
            # convert all iterables to vectors
            return cls.encode_list(value)
        else:
            # Try to convert from numpy if numpy is present
            try:
                import numpy as np

                if isinstance(value, np.number):
                    return str(value)
            except ImportError:
                pass
        raise ValueError(f"Unsupported value for conversion into Julia: {value}")

    @classmethod
    def encode_list(cls, l):
        return "[{}]".format(", ".join(map(cls.encode_value, l)))

    @classmethod
    def encode_items(cls, items):
        def encode_item(item):
            name, value = item
            return f"{cls.encode_value(name)} => {cls.encode_value(value)}"

        return ", ".join(map(encode_item, items))

    @classmethod
    def encode_positional_items(cls, namedlist):
        encoded = ""
        for index, value in enumerate(namedlist):
            encoded += f"{index + 1} => {cls.encode_value(value)}, "
        return encoded

    @classmethod
    def encode_dict(cls, d):
        d = f"Dict({cls.encode_items(d.items())})"
        return d

    @classmethod
    def encode_namedlist(cls, namedlist):
        positional = cls.encode_positional_items(namedlist)
        named = cls.encode_items(namedlist.items())
        source = "Dict("
        if positional:
            source += positional
        if named:
            source += named
        source += ")"
        return source


class BashEncoder:
    """bash docs for associative arrays - https://www.gnu.org/software/bash/manual/html_node/Arrays.html#Arrays"""

    def __init__(
        self,
        namedlists: List[str] = None,
        dicts: List[str] = None,
        prefix: str = "snakemake",
    ):
        """namedlists is a list of strings indicating the snakemake object's member
        variables which are encoded as Namedlist.
        dicts is a list of strings indicating the snakemake object's member variables
        that are encoded as dictionaries.
        Prefix is the prefix for the bash variable name(s) e.g., snakemake_input
        """
        if dicts is None:
            dicts = []
        if namedlists is None:
            namedlists = []
        self.namedlists = namedlists
        self.dicts = dicts
        self.prefix = prefix

    def encode_snakemake(self, smk: Snakemake) -> str:
        """Turn a snakemake object into a collection of bash associative arrays"""
        arrays = []
        main_aa = dict()
        for var in vars(smk):
            val = getattr(smk, var)
            suffix = "params" if var == "_params_store" else var.strip("_").lower()
            if var in self.namedlists:
                aa = f"{self.prefix}_{suffix}={self.encode_namedlist(val)}"
                arrays.append(aa)
            elif var in self.dicts:
                aa = f"{self.prefix}_{suffix}={self.dict_to_aa(val)}"
                arrays.append(aa)
            else:
                main_aa[var] = val

        arrays.append(f"{self.prefix}={self.dict_to_aa(main_aa)}")
        return "\n".join([f"declare -A {aa}" for aa in arrays])

    @staticmethod
    def dict_to_aa(d: dict) -> str:
        """Converts a dictionary to a Bash associative array
        This produces the array component of the variable.
        e.g. ( [var1]=val1 [var2]=val2 )
        to make it a correct bash associative array, you need to name it with
        name=<output of this method>
        """
        s = "( "
        for k, v in d.items():
            # The next replacement is meant for lists, but also gets applied to sub-dicts
            # so that the script will only see a string containing the dict keys.
            # There is no easy way to represent a nested dict in Bash so don't try.
            if isinstance(v, Iterable) and not isinstance(v, str):
                v = " ".join([str(x) for x in v])
            quoted_v = shlex.quote(str(v))
            quoted_k = shlex.quote(str(k))
            s += f"[{quoted_k}]={quoted_v} "

        s += ")"
        return s

    @classmethod
    def encode_namedlist(cls, named_list: iocontainers.Namedlist) -> str:
        """Convert a namedlist into a Bash associative array
        See the comments for dict_to_aa()
        """
        nl_dict = dict()

        # Add the same items keyed by name and also by index
        for i, (name, val) in enumerate(named_list._allitems()):
            if name is not None:
                nl_dict[name] = val
            nl_dict[str(i)] = val

        return cls.dict_to_aa(nl_dict)


class ScriptBase(ABC):
    editable = False

    def __init__(
        self,
        path: SourceFile,
        cache_path: Path,
        source: str,
        is_local: bool,
        run_args: "RunArgs",
        config: Dict,
    ):
        self.path = path
        self.cache_path = cache_path
        self.source = source
        self.is_local = is_local
        self.run_args = run_args
        self.config = config

    def preamble_addendum(self) -> str:
        return ""

    def evaluate(self, edit=False):
        assert not edit or self.editable

        fd = None
        try:
            # generate preamble
            preamble = self.get_preamble()

            # write script
            dir_ = ".snakemake/scripts"
            os.makedirs(dir_, exist_ok=True)

            with tempfile.NamedTemporaryFile(
                suffix="." + self.path.get_filename(), dir=dir_, delete=False
            ) as fd:
                self.write_script(preamble, fd)

            # execute script
            self.execute_script(fd.name, edit=edit)
        except URLError as e:
            raise WorkflowError(e)
        finally:
            if fd and self.run_args.cleanup_scripts:
                os.remove(fd.name)
            else:
                if fd:
                    logger.warning("Not cleaning up %s" % fd.name)
                else:
                    # nothing to clean up (TODO: ??)
                    pass

    @property
    def local_path(self):
        assert self.is_local
        return self.path.get_path_or_uri(secret_free=True)

    @abstractmethod
    def get_preamble(self) -> str: ...

    def write_script(self, preamble, fd):
        fd.write(preamble.encode())
        fd.write(self.source.encode())

    @abstractmethod
    def execute_script(
        self, fname, edit: Optional[NotebookEditMode] = False
    ) -> None: ...

    def _execute_cmd(self, cmd, **kwargs):
        return shell(
            cmd,
            run_args=self.run_args,
            **kwargs,
        )


class PythonScript(ScriptBase):

    @staticmethod
    def _minify_preamble(preamble: str) -> str:
        return textwrap.dedent(preamble).replace("\n", "")

    def preamble_addendum(self) -> str:
        if isinstance(self.path, LocalSourceFile):
            file_override = os.path.realpath(
                self.path.get_path_or_uri(secret_free=True)
            )
        else:
            file_override = self.path.get_path_or_uri(secret_free=True)

        return "__real_file__ = __file__; __file__ = {file_override};".format(
            file_override=repr(file_override)
        )

    def get_preamble(self):
        snakemake = Snakemake(
            input_=self.run_args.input,
            output=self.run_args.output,
            params=self.run_args.params,
            wildcards=self.run_args.wildcards,
            threads=self.run_args.threads,
            resources=self.run_args.resources,
            log=self.run_args.log,
            config=self.config,
            rulename=self.run_args.job_rule.name,
            scriptdir=self.path.get_basedir().get_path_or_uri(secret_free=True),
            bench_iteration=self.run_args.bench_iteration,
        )
        # python 3.14 uses protocol 5 by default, which is not supported in python 3.7 and younger
        # to ensure compatibility, we use protocol 4 here
        snakemake = pickle.dumps(snakemake, protocol=4)
        # Obtain search path for current snakemake module.
        # The module is needed for unpickling in the script.
        # We append it at the end (as a fallback).
        searchpaths = get_snakemake_searchpaths()

        # Add the cache path to the search path so that other cached source files in the same dir
        # can be imported.
        if self.cache_path:
            # TODO handle this in case of container_img, analogously to above
            cache_searchpath = os.path.dirname(self.cache_path)
            if cache_searchpath:
                searchpaths.append(cache_searchpath)
        # For local scripts, add their location to the path in case they use path-based imports
        if self.is_local:
            searchpaths.append(
                self.path.get_basedir().get_path_or_uri(secret_free=True)
            )

        shell_exec = self.run_args.resources.get("shell_exec")
        shell_exec_stmt = (
            ""
            if shell_exec is None
            else f"from snakemake.shell import shell; shell.executable({shell_exec});"
        )

        preamble = f"""
            import sys, pickle;
            sys.path.extend({repr(list(map(str, searchpaths)))});
            snakemake = pickle.loads({snakemake});
            from snakemake.iocontainers import Snakemake;
            is_script = True;
            {shell_exec_stmt}
            {self.preamble_addendum()}
            """
        return "\n".join(
            [
                "######## snakemake preamble start (automatically inserted, do not edit) ########",
                PythonScript._minify_preamble(preamble),
                "######## snakemake preamble end #########\n",
            ]
        )

    def _is_python_env(self) -> bool:
        assert self.run_args.software_env is not None

        is_python_env = self.run_args.software_env.contains_executable("python")
        if ON_WINDOWS and not is_python_env:
            return self.run_args.software_env.contains_executable("python.exe")
        return is_python_env

    def _get_python_version(self):
        # Obtain a clean version string. Using python --version is not reliable, because depending on the distribution
        # stuff may be printed around in unpredictable ways.
        # The code below has to work with python 2.7 as well, therefore it should be written backwards compatible.
        out = self._execute_cmd(
            'python -c "from __future__ import print_function; import sys, json; '
            'print(json.dumps([sys.version_info.major, sys.version_info.minor]))"',
            read=True,
        )
        try:
            return tuple(json.loads(out))
        except ValueError as e:
            raise WorkflowError(
                f"Unable to determine Python version from output '{out}': {e}"
            )

    def execute_script(self, fname, edit: Optional[NotebookEditMode] = False):
        py_exec = sys.executable
        if self.run_args.software_env is not None:
            if self.run_args.software_env.spec.kind == "container":
                # use python from image
                py_exec = "python"
            elif self._is_python_env():
                py_version = self._get_python_version()
                if py_version >= MIN_PY_VERSION:
                    # Python version is new enough, make use of environment
                    # to execute script
                    py_exec = "python"
                else:
                    logger.warning(
                        "Environment defines Python "
                        "version < {0}.{1}. Using Python of the "
                        "main process to execute "
                        "script. Note that this cannot be avoided, "
                        "because the script uses data structures from "
                        "Snakemake which are Python >={0}.{1} "
                        "only.".format(*MIN_PY_VERSION)
                    )

        if ON_WINDOWS:
            # use forward slashes so script command still works even if
            # bash is configured as executable on Windows
            py_exec = py_exec.replace("\\", "/")
        # use the same Python as the running process or the one from the environment
        self._execute_cmd(
            "{py_exec} {fname:q}", py_exec=py_exec, fname=fname, is_python_script=True
        )


class RScript(ScriptBase):
    def get_preamble(self):
        resources = {
            name: value
            for name, value in self.run_args.resources.items()
            if name != "_cores" and name != "_nodes"
        }

        return textwrap.dedent(f"""
        ######## snakemake preamble start (automatically inserted, do not edit) ########
        library(methods)
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
                config = "list",
                rule = "character",
                bench_iteration = "numeric",
                scriptdir = "character",
                source = "function"
            )
        )
        snakemake <- Snakemake(
            input = {REncoder.encode_namedlist(self.run_args.input)},
            output = {REncoder.encode_namedlist(self.run_args.output)},
            params = {REncoder.encode_namedlist(self.run_args.params)},
            wildcards = {REncoder.encode_namedlist(self.run_args.wildcards)},
            threads = {self.run_args.threads},
            log = {REncoder.encode_namedlist(self.run_args.log)},
            resources = {REncoder.encode_namedlist(resources)},
            config = {REncoder.encode_dict(self.config)},
            rule = {REncoder.encode_value(self.run_args.job_rule.name)},
            bench_iteration = {REncoder.encode_numeric(self.run_args.bench_iteration)},
            scriptdir = {REncoder.encode_value(self.path.get_basedir().get_path_or_uri(secret_free=True))},
            source = function(...){{
                old_wd <- getwd()
                on.exit(setwd(old_wd), add = TRUE)

                is_url <- grepl("^https?://", snakemake@scriptdir)
                file <- ifelse(is_url, file.path(snakemake@scriptdir, ...), ...)
                if (!is_url) setwd(snakemake@scriptdir)
                source(file)
            }}
        )

        {self.preamble_addendum()}

        ######## snakemake preamble end #########
        """)

    def execute_script(self, fname, edit: Optional[NotebookEditMode] = False):
        self._execute_cmd("Rscript --vanilla {fname:q}", fname=fname)


class RMarkdown(RScript):

    def write_script(self, preamble, fd):
        # Insert Snakemake object after the RMarkdown header
        code = self.source
        pos = next(itertools.islice(re.finditer(r"---\n", code), 1, 2)).start() + 3
        fd.write(str.encode(code[:pos]))
        preamble = textwrap.dedent("""
            ```{r, echo=FALSE, message=FALSE, warning=FALSE}
            %s
            ```
            """ % preamble)
        fd.write(preamble.encode())
        fd.write(code[pos:].encode())

    def execute_script(self, fname, edit: Optional[NotebookEditMode] = False):
        if len(self.run_args.output) != 1:
            raise WorkflowError(
                "RMarkdown scripts (.Rmd) may only have a single output file."
            )
        out = os.path.abspath(self.run_args.output[0])
        self._execute_cmd(
            'Rscript --vanilla -e \'rmarkdown::render("{fname}", output_file="{out}", quiet=TRUE, knit_root_dir = "{workdir}", params = list(rmd="{fname}"))\'',
            fname=fname,
            out=out,
            workdir=os.getcwd(),
        )


class JuliaScript(ScriptBase):
    def get_preamble(self):
        resources = {
            name: value
            for name, value in self.run_args.resources.items()
            if name != "_cores" and name != "_nodes"
        }
        return textwrap.dedent(f"""
                ######## snakemake preamble start (automatically inserted, do not edit) ########
                struct Snakemake
                    input::Dict
                    output::Dict
                    params::Dict
                    wildcards::Dict
                    threads::Int64
                    log::Dict
                    resources::Dict
                    config::Dict
                    rule::String
                    bench_iteration
                    scriptdir::String
                    #source::Any
                end
                snakemake = Snakemake(
                    {JuliaEncoder.encode_namedlist(self.run_args.input)}, #input::Dict
                    {JuliaEncoder.encode_namedlist(self.run_args.output)}, #output::Dict
                    {JuliaEncoder.encode_namedlist(self.run_args.params)}, #params::Dict
                    {JuliaEncoder.encode_namedlist(self.run_args.wildcards)}, #wildcards::Dict
                    {JuliaEncoder.encode_value(self.run_args.threads)}, #threads::Int64
                    {JuliaEncoder.encode_namedlist(self.run_args.log)}, #log::Dict
                    {JuliaEncoder.encode_namedlist(resources)}, #resources::Dict
                    {JuliaEncoder.encode_dict(self.config)}, #config::Dict
                    {JuliaEncoder.encode_value(self.run_args.job_rule.name)}, #rule::String
                    {JuliaEncoder.encode_value(self.run_args.bench_iteration)}, #bench_iteration::Int64
                    {JuliaEncoder.encode_value(self.path.get_basedir().get_path_or_uri(secret_free=True))}, #scriptdir::String
                    #, #source::Any
                )
                {self.preamble_addendum()}
                ######## snakemake preamble end #########
            """)

    def execute_script(self, fname, edit: Optional[NotebookEditMode] = False):
        self._execute_cmd("julia {fname:q}", fname=fname)


class RustScript(ScriptBase):
    def get_preamble(self):
        # snakemake's namedlists will be encoded as a dict
        # which stores the not-named items at the key "positional"
        # and unpacks named items into the dict
        def encode_namedlist(values):
            values = list(values)
            if len(values) == 0:
                return dict(positional=[])
            positional = [val for key, val in values if not key]
            return dict(
                positional=positional, **{key: val for key, val in values if key}
            )

        snakemake = dict(
            input=encode_namedlist(self.run_args.input._plainstrings()._allitems()),
            output=encode_namedlist(self.run_args.output._plainstrings()._allitems()),
            params=encode_namedlist(self.run_args.params.items()),
            wildcards=encode_namedlist(self.run_args.wildcards.items()),
            threads=self.run_args.threads,
            resources=encode_namedlist(
                {
                    name: value
                    for (name, value) in self.run_args.resources.items()
                    if name != "_cores" and name != "_nodes"
                }.items()
            ),
            log=encode_namedlist(self.run_args.log._plainstrings()._allitems()),
            config=encode_namedlist(self.config.items()),
            rulename=self.run_args.job_rule.name,
            bench_iteration=self.run_args.bench_iteration,
            scriptdir=self.path.get_basedir().get_path_or_uri(secret_free=True),
        )

        json_string = json.dumps(dict(snakemake))

        return textwrap.dedent(f"""
            json_typegen::json_typegen!("Snakemake", r###"{json_string}"###, {{
                "/bench_iteration": {{
                    "use_type": "Option<usize>"
                }},
                "/input/positional": {{
                    "use_type": "Vec<String>"
                }},
                "/output/positional": {{
                    "use_type": "Vec<String>"
                }},
                "/log/positional": {{
                    "use_type": "Vec<String>"
                }},
                "/wildcards/positional": {{
                    "use_type": "Vec<String>"
                }},
            }});

            pub struct Iter<'a, T>(std::slice::Iter<'a, T>);
            impl<'a, T> Iterator for Iter<'a, T> {{
                type Item = &'a T;

                fn next(&mut self) -> Option<Self::Item> {{
                    self.0.next()
                }}
            }}
            macro_rules! impl_iter {{
                ($($s:ty),+) => {{
                    $(
                        impl IntoIterator for $s {{
                            type Item = String;
                            type IntoIter = std::vec::IntoIter<Self::Item>;

                            fn into_iter(self) -> Self::IntoIter {{
                                self.positional.into_iter()
                            }}
                        }}

                        impl<'a> IntoIterator for &'a $s {{
                            type Item = &'a String;
                            type IntoIter = Iter<'a, String>;

                            fn into_iter(self) -> Self::IntoIter {{
                                Iter(self.positional.as_slice().into_iter())
                            }}
                        }}
                    )+
                }};
            }}

            macro_rules! impl_index {{
                ($($s:ty),+) => {{
                    $(
                    impl std::ops::Index<usize> for $s {{
                        type Output = String;

                        fn index(&self, index: usize) -> &Self::Output {{
                            &self.positional[index]
                        }}
                    }}
                    )+
                }}
            }}


            impl_iter!(Input, Output, Wildcards, Log);
            impl_index!(Input, Output, Wildcards, Log);

            impl Snakemake {{
                #[allow(dead_code)]
                fn redirect_stderr<P: AsRef<std::path::Path>>(
                    &self,
                    path: P,
                ) -> anyhow::Result<gag::Redirect<std::fs::File>> {{
                    let log = std::fs::OpenOptions::new()
                        .truncate(true)
                        .read(true)
                        .create(true)
                        .write(true)
                        .open(path)?;
                    Ok(gag::Redirect::stderr(log)?)
                }}

                #[allow(dead_code)]
                fn redirect_stdout<P: AsRef<std::path::Path>>(
                    &self,
                    path: P,
                ) -> anyhow::Result<gag::Redirect<std::fs::File>> {{
                    let log = std::fs::OpenOptions::new()
                        .truncate(true)
                        .read(true)
                        .create(true)
                        .write(true)
                        .open(path)?;
                    Ok(gag::Redirect::stdout(log)?)
                }}
            }}

            lazy_static::lazy_static! {{
                // https://github.com/rust-lang-nursery/lazy-static.rs/issues/153
                #[allow(non_upper_case_globals)]
                static ref snakemake: Snakemake = {{
                    let s: Snakemake = serde_json::from_str(r###"{json_string}"###).expect("Failed parsing snakemake JSON");
                    s
                }};
            }}
            {self.preamble_addendum()}
            """)

    def write_script(self, preamble, fd):
        content = self.combine_preamble_and_source(preamble)
        fd.write(content.encode())

    def execute_script(self, fname, edit: Optional[NotebookEditMode] = False):
        deps = self.default_dependencies()
        ftrs = self.default_features()
        self._execute_cmd(
            "rust-script -d {deps} {fname:q} -- --features {ftrs}",
            fname=fname,
            deps=deps,
            ftrs=ftrs,
        )

    def combine_preamble_and_source(self, preamble: str) -> str:
        """The manifest info needs to be moved to before the preamble.
        Also, because rust-scipt relies on inner docs, there can't be an empty line
        between the manifest and preamble.
        """
        manifest, src = RustScript.extract_manifest(self.source)
        return manifest + preamble.lstrip("\r\n") + src

    @staticmethod
    def default_dependencies() -> str:
        return " -d ".join(
            [
                "anyhow=1",
                "serde_json=1.0",
                "serde=1.0",
                "serde_derive=1.0",
                "lazy_static=1.4",
                "json_typegen=0.6",
                "gag=1",
            ]
        )

    @staticmethod
    def default_features() -> str:
        return ",".join(["serde/derive"])

    @staticmethod
    def extract_manifest(source: str) -> Tuple[str, str]:
        # we have no need for the shebang for now given the way we run the script
        _, src = RustScript._strip_shebang(source)
        manifest, src = RustScript._strip_manifest(src)

        return manifest, src

    @staticmethod
    def _strip_shebang(src: str) -> Tuple[str, str]:
        """From https://github.com/fornwall/rust-script/blob/ce508bad02a11d574657d2f1debf7e73fca2bf6e/src/manifest.rs#L312-L320"""
        rgx = re.compile(r"^#![^\[].*?(\r\n|\n)")
        return strip_re(rgx, src)

    @staticmethod
    def _strip_manifest(src: str) -> Tuple[str, str]:
        """From https://github.com/fornwall/rust-script/blob/ce508bad02a11d574657d2f1debf7e73fca2bf6e/src/manifest.rs#L405-L411"""
        manifest, remainder = RustScript._strip_single_line_manifest(src)
        if not manifest:
            manifest, remainder = RustScript._strip_code_block_manifest(src)
        return manifest, remainder

    @staticmethod
    def _strip_single_line_manifest(src: str) -> Tuple[str, str]:
        """From https://github.com/fornwall/rust-script/blob/ce508bad02a11d574657d2f1debf7e73fca2bf6e/src/manifest.rs#L618-L632"""
        rgx = re.compile(r"^\s*//\s*cargo-deps\s*:(.*?)(\r\n|\n)", flags=re.IGNORECASE)
        return strip_re(rgx, src)

    @staticmethod
    def _strip_code_block_manifest(src: str) -> Tuple[str, str]:
        """From https://github.com/fornwall/rust-script/blob/ce508bad02a11d574657d2f1debf7e73fca2bf6e/src/manifest.rs#L634-L664
        We need to find the first `/*!` or `//!` that *isn't* preceded by something
        that would make it apply to anything other than the create itself. Because we
        can't do this accurately, we'll just require that the doc comment is the
        *first* thing in the file (after the optional shebang, which should already
        have been stripped).
        """
        crate_comment_re = re.compile(
            r"^\s*(/\*!|//([!/]))(.*?)(\r\n|\n)", flags=re.MULTILINE
        )
        # does src start with a create comment?
        match = crate_comment_re.match(src)
        if not match:
            return "", src
        end_of_comment = match.end()
        # find end of create comment
        while match is not None:
            end_of_comment = match.end()
            match = crate_comment_re.match(src, pos=end_of_comment)

        crate_comment = src[:end_of_comment]
        found_code_block_open = False
        code_block_open_re = re.compile(r"```\s*cargo")
        found_code_block_close = False
        code_block_close_re = re.compile(r"```")
        for line in crate_comment.splitlines():
            if not found_code_block_open:
                m = code_block_open_re.search(line)
                if m:
                    found_code_block_open = True
            else:
                m = code_block_close_re.search(line)
                if m:
                    found_code_block_close = True
                    break

        crate_comment_has_manifest = found_code_block_open and found_code_block_close
        if crate_comment_has_manifest:
            return crate_comment, src[end_of_comment:]
        else:
            return "", src


class BashScript(ScriptBase):
    def get_preamble(self):
        snakemake = Snakemake(
            input_=self.run_args.input,
            output=self.run_args.output,
            params=self.run_args.params,
            wildcards=self.run_args.wildcards,
            threads=self.run_args.threads,
            resources=self.run_args.resources,
            log=self.run_args.log,
            config=self.config,
            rulename=self.run_args.job_rule.name,
            scriptdir=self.path.get_basedir().get_path_or_uri(secret_free=True),
            bench_iteration=self.run_args.bench_iteration,
        )

        namedlists = [
            "input",
            "output",
            "log",
            "resources",
            "wildcards",
            "_params_store",
        ]
        dicts = ["config"]
        encoder = BashEncoder(namedlists=namedlists, dicts=dicts)
        preamble = encoder.encode_snakemake(snakemake)
        return f"{preamble}\n{self.preamble_addendum()}"

    def write_script(self, preamble, fd):
        content = self.combine_preamble_and_source(preamble)
        fd.write(content.encode())

    def combine_preamble_and_source(self, preamble: str):
        rgx = re.compile(r"^#![^\[].*?(\r\n|\n)")
        shebang, source = strip_re(rgx, self.source)
        if not shebang:
            shebang = r"#!/usr/bin/env bash"

        return "\n".join([shebang, preamble, source])

    def execute_script(self, fname, edit: Optional[NotebookEditMode] = False):
        self._execute_cmd("bash {fname:q}", fname=fname)


class XonshScript(PythonScript):
    def execute_script(self, fname, edit: Optional[NotebookEditMode] = False):
        self._execute_cmd(
            "xonsh -DRAISE_SUBPROC_ERROR=true -DXONSH_SHOW_TRACEBACK=true {fname:q}",
            fname=fname,
        )


class HyScript(PythonScript):
    def write_script(self, preamble, fd):
        fd.write(f"(pys #[[{preamble}]])".encode())
        fd.write(self.source.encode())

    def execute_script(self, fname, edit: Optional[NotebookEditMode] = False):
        self._execute_cmd("hy {fname:q}", fname=fname)


def strip_re(regex: Pattern, s: str) -> Tuple[str, str]:
    """Strip a substring matching a regex from a string and return the stripped part
    and the remainder of the original string.
    Returns an empty string and the original string if the regex is not found
    """
    rgx = re.compile(regex)
    match = rgx.search(s)
    if match:
        head, tail = s[: match.end()], s[match.end() :]
    else:
        head, tail = "", s

    return head, tail


def get_source(
    path,
    sourcecache: SourceCache,
    basedir=None,
    wildcards=None,
    params=None,
):
    if wildcards is not None and params is not None:
        if isinstance(path, SourceFile):
            path = path.format(wildcards=wildcards, params=params)
        else:
            path = infer_source_file(format(path, wildcards=wildcards, params=params))

    if basedir is not None:
        basedir = infer_source_file(basedir)

    source_file = infer_source_file(path, basedir)
    with sourcecache.open(source_file) as f:
        source = f.read()

    language = get_language(source_file, source)

    is_local = isinstance(source_file, LocalSourceFile)

    return source_file, source, language, is_local, sourcecache.get_path(source_file)


def get_language(source_file, source):
    import nbformat

    filename = source_file.get_filename()

    language = None
    if filename.endswith(".py"):
        language = "python"
    elif filename.endswith(".ipynb"):
        language = "jupyter"
    elif filename.endswith(".R"):
        language = "r"
    elif filename.endswith(".Rmd"):
        language = "rmarkdown"
    elif filename.endswith(".jl"):
        language = "julia"
    elif filename.endswith(".rs"):
        language = "rust"
    elif filename.endswith(".sh"):
        language = "bash"
    elif filename.endswith(".xsh"):
        language = "xonsh"
    elif filename.endswith(".hy"):
        language = "hy"

    # detect kernel language for Jupyter Notebooks
    if language == "jupyter":
        nb = nbformat.reads(source, as_version=nbformat.NO_CONVERT)
        try:
            kernel_language = nb["metadata"]["language_info"]["name"]
        except KeyError as e:
            raise WorkflowError(
                "Notebook metadata is corrupt. Please delete notebook "
                "and recreate it via --edit-notebook."
            )

        language += "_" + kernel_language.lower()

    return language


def script(
    path,
    run_args: "RunArgs",
    config: Dict,
):
    """
    Load a script from the given basedir + path and execute it.
    """
    if isinstance(path, Path):
        path = str(path)

    path, source, language, is_local, cache_path = get_source(
        path,
        SourceCache(run_args.cache_path, run_args.runtime_cache_path),
        run_args.basedir,
        run_args.wildcards,
        run_args.params,
    )

    exec_class = {
        "python": PythonScript,
        "r": RScript,
        "rmarkdown": RMarkdown,
        "julia": JuliaScript,
        "rust": RustScript,
        "bash": BashScript,
        "xonsh": XonshScript,
        "hy": HyScript,
    }.get(language, None)
    if exec_class is None:
        raise ValueError(
            "Script must be one of the following filetypes: [.py .R .Rmd .jl .rs .sh .xsh .hy]"
        )

    executor = exec_class(
        path=path,
        cache_path=cache_path,
        source=source,
        is_local=is_local,
        run_args=run_args,
        config=config,
    )
    executor.evaluate()
