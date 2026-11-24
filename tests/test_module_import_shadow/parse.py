# A module whose name collides with a snakemake internal (workflow.py imports
# snakemake.parser.parse). Importing it inside a Snakefile must not break the
# 'module' keyword or any other workflow statement.
VALUE = 1
