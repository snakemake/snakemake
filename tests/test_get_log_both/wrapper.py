from snakemake.shell import shell
shell.executable("bash")

log = snakemake.log_fmt_shell(append=True)
shell(
    """
      cat {snakemake.input} > {snakemake.output}
      (>&2 echo "a stderr message") {log}
      (echo "a stdout message") {log}
      """
)
