rule all:
  input: "python-version.txt"

rule py_emit_version:
  output: "python-version.txt"
  conda: "envs/py39.yaml"
  shell:
    """
    python -VV > {output}
    """
