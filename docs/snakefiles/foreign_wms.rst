
.. _snakefiles-foreign-wms:

===============================================
Integrating foreign workflow management systems
===============================================

Snakemake 6.2 and later allows to hand over execution steps to other workflow management systems.
By this, it is possible to make use of workflows written for other systems, while performing any further pre- or postprocessing within Snakemake.
Such a handover is indicated with the ``handover`` directive.
Consider the following example:

.. code-block:: python

    rule chipseq_pipeline:
        input:
            input="design.csv",
            fasta="data/genome.fasta",
            gtf="data/genome.gtf",
        output:
            "multiqc/broadPeaks/multiqc_report.html",
        params:
            pipeline="nf-core/chipseq",
            revision="1.2.1",
            profile=["conda"],
        handover: True
        wrapper:
            "0.74.0/utils/nextflow"

Here, the workflow is executed as usual until this rule is reached.
Then, Snakemake passes all resources to the nextflow workflow management system, which generates certain files.
The rule is executed as a :ref:`local rule <snakefiles-local-rules>`, meaning that it would not be submitted to a cluster or cloud system by Snakemake.
Instead, the invoked other workflow management system is responsible for that.
E.g., in case of `Nextflow <https://nextflow.io>`_, submission behavior can be configured via a ``nextflow.conf`` file or environment variables.
After the step is done, Snakemake continues execution with the output files produced by the foreign workflow.
