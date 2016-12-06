.. _user_manual-writing_snakefiles:

==================
Writing Snakefiles
==================

In Snakemake, workflows are specified as Snakefiles. Inspired by GNU Make, a Snakefile contains rules, that denote how to create output files from input files.
Dependencies between rules are handled implicitly, by matching filenames of input files against output files. Thereby wildcards can be used to write general rules.
