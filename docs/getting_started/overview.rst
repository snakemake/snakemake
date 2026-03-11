This tutorial introduces the text-based workflow system Snakemake.

Snakemake follows the GNU Make paradigm: workflows are defined in terms of rules that define how to create output files from input files.
Dependencies between the rules are determined automatically, creating a DAG (directed acyclic graph) of jobs that can be automatically parallelized.

Snakemake sets itself apart from other text-based workflow systems in the following way.
Hooking into the Python interpreter, Snakemake offers a definition language that is an extension of Python with syntax to define rules and workflow specific properties.
This allows Snakemake to combine the flexibility of a plain scripting language with a pythonic workflow definition.
The Python language is known to be concise yet readable and can appear almost like pseudo-code.
The syntactic extensions provided by Snakemake maintain this property for the definition of the workflow.
Further, Snakemake’s scheduling algorithm can be constrained by priorities, provided cores and customizable resources and it provides a generic support for distributed computing (e.g., cluster or batch systems).
Hence, a Snakemake workflow scales without modification from single core workstations and multi-core servers to cluster or batch systems.
Finally, Snakemake integrates with the package manager Conda and the container engine Singularity such that defining the software stack becomes part of the workflow itself.

Historically, Snakemake has been developed with the domain of Bioinformatics in mind.
However, Snakemake is a general-purpose workflow management system for any discipline.


