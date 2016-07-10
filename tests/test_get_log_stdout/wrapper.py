from snakemake.shell import shell
log = snakemake.get_log(stderr=False, append=True)
shell('''
      cat {snakemake.input} > {snakemake.output}
      (>&2 echo "a stderr message") {log}
      (echo "a stdout message") {log}
      ''')
