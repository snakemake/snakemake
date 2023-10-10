from snakemake.shell import shell
shell.executable("bash")

import sys

sys.stderr.write("a stderr message from the wrapper\n")
sys.stdout.write("a stdout message from the wrapper\n")

print("A message goes to no log file")

log = snakemake.log_fmt_shell(append=False)
shell(
    """
      cat {snakemake.input} > {snakemake.output}
      (>&2 echo "a stderr message") {log}
      (echo "a stdout message") {log}
      """
)
