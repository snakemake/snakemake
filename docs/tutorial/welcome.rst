.. _tutorial-welcome:

Snakemake Tutorial
==================

.. _Snakemake: http://snakemake.bitbucket.org
.. _Snakemake homepage: http://snakemake.bitbucket.org
.. _GNU Make: https://www.gnu.org/software/make
.. _Python: http://www.python.org
.. _BWA: http://bio-bwa.sourceforge.net
.. _SAMtools: http://www.htslib.org
.. _BCFtools: http://www.htslib.org
.. _Pandas: http://pandas.pydata.org
.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _Conda: http://conda.pydata.org
.. _Bash: http://www.tldp.org/LDP/Bash-Beginners-Guide/html
.. _Atom: https://atom.io
.. _Anaconda: https://anaconda.org
.. _Graphviz: http://www.graphviz.org
.. _RestructuredText: http://docutils.sourceforge.net/rst.html
.. _data URI: https://developer.mozilla.org/en-US/docs/Web/HTTP/data_URIs
.. _Documentation: https://bitbucket.org/snakemake/snakemake/wiki/Documentation
.. _JSON: http://json.org
.. _YAML: http://yaml.org
.. _DRMAA: http://www.drmaa.org
.. _FAQ: https://bitbucket.org/snakemake/snakemake/wiki/FAQ
.. _rpy2: http://rpy.sourceforge.net
.. _R: https://www.r-project.org
.. _Rscript: https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html
.. _cluster configuration: https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-cluster-configuration
.. _script section in the Documentation: https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-external-scripts
.. _PyYAML: http://pyyaml.org
.. _Docutils: http://docutils.sourceforge.net
.. _Bioconda: https://bioconda.github.io
.. _Vagrant: https://www.vagrantup.com
.. _Vagrant Documentation: https://docs.vagrantup.com
.. _Blogpost: http://blog.osteel.me/posts/2015/01/25/how-to-use-vagrant-on-windows.html
.. _slides: http://slides.com/johanneskoester/deck-1

This tutorial introduces the text-based workflow system Snakemake_.
Snakemake follows the `GNU Make`_ paradigm: workflows are defined in terms of rules that define how to create output files from input files.
Dependencies between the rules are determined automatically, creating a DAG (directed acyclic graph) of jobs that can be automatically parallelized.

Snakemake sets itself apart from existing text-based workflow systems in the following way.
Hooking into the Python interpreter, Snakemake offers a definition language that is an extension of Python_ with syntax to define rules and workflow specific properties.
This allows to combine the flexibility of a plain scripting language with a pythonic workflow definition.
The Python language is known to be concise yet readable and can appear almost like pseudo-code.
The syntactic extensions provided by Snakemake maintain this property for the definition of the workflow.
Further, Snakemakes scheduling algorithm can be constrained by priorities, provided cores and customizable resources and it provides a generic support for distributed computing (e.g., cluster or batch systems).
Hence, a Snakemake workflow scales without modification from single core workstations and multi-core servers to cluster or batch systems.

While the examples presented here come from Bioinformatics, Snakemake is considered a general-purpose workflow management system for any discipline.

Also have a look at the corresponding slides_.
