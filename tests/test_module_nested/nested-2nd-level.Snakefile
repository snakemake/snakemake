from snakemake.utils import min_version
min_version("6.0")

module module_3:
    snakefile: "nested-3rd-level.Snakefile"
    prefix: "nested-3rd-level"

rule default:
    input:
       # rules.nested_3rd_level_default.input
       "done"
    default_target: True

rule do:
    input:
      "nested-3rd-level/done"
    output:
      "done"
    shell:
      "touch {output}"

use rule * from module_3 exclude default as nested_3rd_level_*


