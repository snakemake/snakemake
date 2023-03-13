from snakemake.utils import min_version
min_version("6.0")

module module_2:
    snakefile: "nested-2nd-level.Snakefile"
    prefix: "nested-2nd-level"

rule do:
    input:
      "nested-2nd-level/done"
    output:
      "done"
    shell:
      "touch {output}"

rule default:
    input:
       # rules.nested_2nd_level_default.input
       "done"
    default_target: True

use rule * from module_2 exclude default as nested_2nd_level_*
