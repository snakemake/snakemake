__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2020, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

python_script_template = """
        ######## snakemake preamble start (automatically inserted, do not edit) ########
        import sys, os; sys.path.extend([{searchpath}]); sys.path.insert(0, {sourcedir});
        class Snakemake:
            def __init__(self, input, output, params, wildcards, threads, log, resources, config, rule, bench_iteration, scriptdir, source=None):
                self.set_listvar("input", input)
                self.set_listvar("params", params)
                self.set_listvar("output", output)
                self.wildcards = wildcards
                self.threads = threads
                self.log = log
                self.resources = resources
                self.config = config
                self.rule = rule
                self.bench_iteration = bench_iteration
                self.scriptdir = scriptdir
                if source is not None:
                    self.source = source

            def set_listvar(self, name, items):
                class Item:
                    def __init__(self, entries):
                        self.__dict__.update(entries)

                    def get(self, key, default=None):
                        if hasattr(self, key):
                            return getattr(self, key)
                        return default

                    def __getitem__(self, key, default=None):
                        try:
                            return self.__dict__[key]
                        except:
                            return default

                if isinstance(items, dict):
                    structure = Item(items)
                    setattr(self, name, structure)
                else:
                    setattr(self, name, items) 

            def source(self):
                wd = os.getcwd()
                os.chdir(self.scriptdir)
                sys.path.insert(self.scriptdir)
                os.chdir(wd)
                
            def log_fmt_shell(self, stdout=True, stderr=True, append=False):
                if not self.log:
                    return ""
                if stdout == True and stderr == True and append == True:
                    return " >> {{0}} 2>&1".format(self.log[0])
                elif stdout == True and stderr == False and append == True:
                    return " >> {{0}}".format(self.log[0])
                elif stdout == False and stderr == True and append == True:
                    return " 2>> {{0}}".format(self.log[0])
                elif stdout == True and stderr == True and append == False:
                    return " > {{0}} 2>&1".format(self.log[0])
                elif stdout == True and stderr == False and append == False:
                    return " > {{0}}".format(self.log[0])
                else: # False True False
                    return " 2> {{0}}".format(self.log[0])

        snakemake = Snakemake(
            input = {input_},
            output = {output},
            params = {params},
            wildcards = {wildcards},
            threads = {threads},
            log = {log},
            resources = {resources},
            config = {config},
            rule = {rule},
            bench_iteration = {bench_iteration},
            scriptdir = {scriptdir},
        )
        def shell(cmd):
            import inspect, os
            context = inspect.currentframe().f_back.f_locals
            os.system(cmd.format(**context))
        {preamble_addendum}

        ######## snakemake preamble end #########
        """

r_script_template = """
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
            input = {},
            output = {},
            params = {},
            wildcards = {},
            threads = {},
            log = {},
            resources = {},
            config = {},
            rule = {},
            bench_iteration = {},
            scriptdir = {},
            source = function(...){{
                wd <- getwd()
                setwd(snakemake@scriptdir)
                source(...)
                setwd(wd)
            }}
        )
        {preamble_addendum}

        ######## snakemake preamble end #########
        """

rmarkdown_script_template = """
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
            input = {},
            output = {},
            params = {},
            wildcards = {},
            threads = {},
            log = {},
            resources = {},
            config = {},
            rule = {},
            bench_iteration = {},
            scriptdir = {},
            source = function(...){{
                wd <- getwd()
                setwd(snakemake@scriptdir)
                source(...)
                setwd(wd)
            }}
        )

        ######## snakemake preamble end #########
        """

julia_script_template = """
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
                    {}, #input::Dict
                    {}, #output::Dict
                    {}, #params::Dict
                    {}, #wildcards::Dict
                    {}, #threads::Int64
                    {}, #log::Dict
                    {}, #resources::Dict
                    {}, #config::Dict
                    {}, #rule::String
                    {}, #bench_iteration::Int64
                    {}, #scriptdir::String
                    #, #source::Any
                )
                ######## snakemake preamble end #########
                """
