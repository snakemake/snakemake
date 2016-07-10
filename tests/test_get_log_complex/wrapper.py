from snakemake.shell import shell
initial_log = snakemake.get_log()
stdout_log = snakemake.get_log(stderr=False, append=True)
stderr_log = snakemake.get_log(stdout=False, append=True)
shell('''
      cat {snakemake.input} > {snakemake.output}
      echo "should not appear since next line truncates" {initial_log}
      echo "first line" {initial_log}
      (>&2 echo "a stderr message") {stderr_log}
      (echo "a stdout message") {stdout_log}
      ''')
