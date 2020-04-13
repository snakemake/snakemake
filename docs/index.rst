.. _manual-main:

=========
Snakemake
=========

.. image:: https://img.shields.io/conda/dn/bioconda/snakemake.svg?label=Bioconda
    :target: https://bioconda.github.io/recipes/snakemake/README.html

.. image:: https://img.shields.io/pypi/pyversions/snakemake.svg
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/snakemake.svg
    :target: https://pypi.python.org/pypi/snakemake

.. image:: https://img.shields.io/docker/cloud/build/snakemake/snakemake
       :target: https://hub.docker.com/r/snakemake/snakemake

.. image:: https://github.com/snakemake/snakemake/workflows/CI/badge.svg?branch=master
    :target: https://github.com/snakemake/snakemake/actions?query=branch%3Amaster+workflow%3ACI

.. image:: https://img.shields.io/badge/stack-overflow-orange.svg
    :target: https://stackoverflow.com/questions/tagged/snakemake

.. image:: https://img.shields.io/twitter/follow/johanneskoester.svg?style=social&label=Follow
    :target: https://twitter.com/search?l=&q=%23snakemake%20from%3Ajohanneskoester

.. image:: https://img.shields.io/github/stars/snakemake/snakemake?style=social
    :alt: GitHub stars
    :target: https://github.com/snakemake/snakemake/stargazers

.. .. raw:: html
          <span class="__dimensions_badge_embed__" data-doi="https://doi.org/10.1093/bioinformatics/bts480" data-legend="always" data-style="large_rectangle"></span><script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>

The Snakemake workflow management system is a tool to create **reproducible and scalable** data analyses.
Workflows are described via a human readable, Python based language.
They can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition.
Finally, Snakemake workflows can entail a description of required software, which will be automatically deployed to any execution environment.

Snakemake is **highly popular** with, `~3 new citations per week <https://badge.dimensions.ai/details/id/pub.1018944052>`_.


.. _manual-quick_example:

-------------
Quick Example
-------------

Snakemake workflows are essentially Python scripts extended by declarative code to define **rules**.
Rules describe how to create **output files** from **input files**.

.. code-block:: python

    rule targets:
        input:
            "plots/myplot.pdf"

    rule transform:
        input:
            "raw/{dataset}.csv"
        output:
            "transformed/{dataset}.csv"
        singularity:
            "docker://somecontainer:v1.0"
        shell:
            "somecommand {input} {output}"

    rule aggregate_and_plot:
        input:
            expand("transformed/{dataset}.csv", dataset=[1, 2])
        output:
            "plots/myplot.pdf"
        conda:
            "envs/matplotlib.yaml"
        script:
            "scripts/plot.py"


* Similar to GNU Make, you specify targets in terms of a pseudo-rule at the top.
* For each target and intermediate file, you create rules that define how they are created from input files.
* Snakemake determines the rule dependencies by matching file names.
* Input and output files can contain multiple named wildcards.
* Rules can either use shell commands, plain Python code or external Python or R scripts to create output files from input files.
* Snakemake workflows can be easily executed on **workstations**, **clusters**, **the grid**, and **in the cloud** without modification. The job scheduling can be constrained by arbitrary resources like e.g. available CPU cores, memory or GPUs.
* Snakemake can automatically deploy required software dependencies of a workflow using `Conda <https://conda.io>`_ or `Singularity <https://sylabs.io/docs/>`_.
* Snakemake can use Amazon S3, Google Storage, Dropbox, FTP, WebDAV, SFTP and iRODS to access input or output files and further access input files via HTTP and HTTPS.


.. _main-getting-started:

---------------
Getting started
---------------

To get a first impression, see our `introductory slides <https://slides.com/johanneskoester/snakemake-short>`_ or watch the `live demo video <https://youtu.be/hPrXcUUp70Y>`_.
News about Snakemake are published via `Twitter <https://twitter.com/search?l=&q=%23snakemake%20from%3Ajohanneskoester>`_.
To learn Snakemake, please do the :ref:`tutorial`, and see the :ref:`FAQ <project_info-faq>`.

.. _main-support:

-------
Support
-------

* For releases, see :ref:`Changelog <changelog>`.
* Check :ref:`frequently asked questions (FAQ) <project_info-faq>`.
* In case of questions, please post on `stack overflow <https://stackoverflow.com/questions/tagged/snakemake>`_.
* To discuss with other Snakemake users, you can use the `mailing list <https://groups.google.com/forum/#!forum/snakemake>`_. **Please do not post questions there. Use stack overflow for questions.**
* For bugs and feature requests, please use the `issue tracker <https://github.com/snakemake/snakemake/issues>`_.
* For contributions, visit Snakemake on `Github <https://github.com/snakemake/snakemake>`_ and read the :ref:`guidelines <project_info-contributing>`.

--------
Citation
--------

`Köster, Johannes and Rahmann, Sven. "Snakemake - A scalable bioinformatics workflow engine". Bioinformatics 2012. <https://bioinformatics.oxfordjournals.org/content/28/19/2520>`_

See :doc:`Citations <project_info/citations>` for more information.

---------
Resources
---------

`Snakemake Wrappers Repository <https://snakemake-wrappers.readthedocs.org>`_
    The Snakemake Wrapper Repository is a collection of reusable wrappers that allow to quickly use popular tools from Snakemake rules and workflows.

`Snakemake Workflows Project <https://github.com/snakemake-workflows/docs>`_
    This project provides a collection of high quality modularized and re-usable workflows.
    The provided code should also serve as a best-practices of how to build production ready workflows with Snakemake.
    Everybody is invited to contribute.

`Snakemake Profiles Project <https://github.com/snakemake-profiles/doc>`_
    This project provides Snakemake configuration profiles for various execution environments.
    Please consider contributing your own if it is still missing.

`Bioconda <https://bioconda.github.io/>`_
    Bioconda can be used from Snakemake for creating completely reproducible workflows by defining the used software versions and providing binaries.


.. project_info-publications_using:

----------------------------
Publications using Snakemake
----------------------------

In the following you find an **incomplete list** of publications making use of Snakemake for their analyses.
Please consider to add your own.

* Kuzniar et al. 2020. `sv-callers: a highly portable parallel workflow for structural variant detection in whole-genome sequence data <https://doi.org/10.7717/peerj.8214>`_. PeerJ.
* Doris et al. 2018. `Spt6 is required for the fidelity of promoter selection <https://doi.org/10.1016/j.molcel.2018.09.005>`_. Molecular Cell.
* Karlsson et al. 2018. `Four evolutionary trajectories underlie genetic intratumoral variation in childhood cancer <https://www.nature.com/articles/s41588-018-0131-y>`_. Nature Genetics.
* Planchard et al. 2018. `The translational landscape of Arabidopsis mitochondria <https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gky489/5033161>`_. Nucleic acids research.
* Schult et al. 2018. `Effect of UV irradiation on Sulfolobus acidocaldarius and involvement of the general transcription factor TFB3 in the early UV response <https://academic.oup.com/nar/article/46/14/7179/5047281>`_. Nucleic acids research.
* Goormaghtigh et al. 2018. `Reassessing the Role of Type II Toxin-Antitoxin Systems in Formation of Escherichia coli Type II Persister Cells <https://mbio.asm.org/content/mbio/9/3/e00640-18.full.pdf>`_. mBio.
* Ramirez et al. 2018. `Detecting macroecological patterns in bacterial communities across independent studies of global soils <https://www.nature.com/articles/s41564-017-0062-x>`_. Nature microbiology.
* Amato et al. 2018. `Evolutionary trends in host physiology outweigh dietary niche in structuring primate gut microbiomes <https://www.nature.com/articles/s41396-018-0175-0>`_. The ISME journal.
* Uhlitz et al. 2017. `An immediate–late gene expression module decodes ERK signal duration <https://msb.embopress.org/content/13/5/928>`_. Molecular Systems Biology.
* Akkouche et al. 2017. `Piwi Is Required during Drosophila Embryogenesis to License Dual-Strand piRNA Clusters for Transposon Repression in Adult Ovaries <https://www.sciencedirect.com/science/article/pii/S1097276517302071>`_. Molecular Cell.
* Beatty et al. 2017. `Giardia duodenalis induces pathogenic dysbiosis of human intestinal microbiota biofilms <https://www.ncbi.nlm.nih.gov/pubmed/28237889>`_. International Journal for Parasitology.
* Meyer et al. 2017. `Differential Gene Expression in the Human Brain Is Associated with Conserved, but Not Accelerated, Noncoding Sequences <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5400397/>`_. Molecular Biology and Evolution.
* Lonardo et al. 2017. `Priming of soil organic matter: Chemical structure of added compounds is more important than the energy content <https://www.sciencedirect.com/science/article/pii/S0038071716304539>`_. Soil Biology and Biochemistry.
* Beisser et al. 2017. `Comprehensive transcriptome analysis provides new insights into nutritional strategies and phylogenetic relationships of chrysophytes <https://peerj.com/articles/2832/>`_. PeerJ.
* Piro et al 2017. `MetaMeta: integrating metagenome analysis tools to improve taxonomic profiling <https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0318-y>`_. Microbiome.
* Dimitrov et al 2017. `Successive DNA extractions improve characterization of soil microbial communities <https://peerj.com/articles/2915/>`_. PeerJ.
* de Bourcy et al. 2016. `Phylogenetic analysis of the human antibody repertoire reveals quantitative signatures of immune senescence and aging <https://www.pnas.org/content/114/5/1105.short>`_. PNAS.
* Bray et al. 2016. `Near-optimal probabilistic RNA-seq quantification <https://www.nature.com/nbt/journal/v34/n5/abs/nbt.3519.html>`_. Nature Biotechnology.
* Etournay et al. 2016. `TissueMiner: a multiscale analysis toolkit to quantify how cellular processes create tissue dynamics <https://elifesciences.org/content/5/e14334>`_. eLife Sciences.
* Townsend et al. 2016. `The Public Repository of Xenografts Enables Discovery and Randomized Phase II-like Trials in Mice <https://www.cell.com/cancer-cell/abstract/S1535-6108%2816%2930090-3>`_. Cancer Cell.
* Burrows et al. 2016. `Genetic Variation, Not Cell Type of Origin, Underlies the Majority of Identifiable Regulatory Differences in iPSCs <https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005793>`_. PLOS Genetics.
* Ziller et al. 2015. `Coverage recommendations for methylation analysis by whole-genome bisulfite sequencing <https://www.nature.com/nmeth/journal/v12/n3/full/nmeth.3152.html>`_. Nature Methods.
* Li et al. 2015. `Quality control, modeling, and visualization of CRISPR screens with MAGeCK-VISPR <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0843-6>`_. Genome Biology.
* Schmied et al. 2015. `An automated workflow for parallel processing of large multiview SPIM recordings <https://bioinformatics.oxfordjournals.org/content/32/7/1112>`_. Bioinformatics.
* Chung et al. 2015. `Whole-Genome Sequencing and Integrative Genomic Analysis Approach on Two 22q11.2 Deletion Syndrome Family Trios for Genotype to Phenotype Correlations <https://onlinelibrary.wiley.com/doi/10.1002/humu.22814/full>`_. Human Mutation.
* Kim et al. 2015. `TUT7 controls the fate of precursor microRNAs by using three different uridylation mechanisms <https://emboj.embopress.org/content/34/13/1801.long>`_. The EMBO Journal.
* Park et al. 2015. `Ebola Virus Epidemiology, Transmission, and Evolution during Seven Months in Sierra Leone <https://doi.org/10.1016/j.cell.2015.06.007>`_. Cell.
* Břinda et al. 2015. `RNF: a general framework to evaluate NGS read mappers <https://bioinformatics.oxfordjournals.org/content/early/2015/09/30/bioinformatics.btv524>`_. Bioinformatics.
* Břinda et al. 2015. `Spaced seeds improve k-mer-based metagenomic classification <https://bioinformatics.oxfordjournals.org/content/early/2015/08/10/bioinformatics.btv419>`_. Bioinformatics.
* Spjuth et al. 2015. `Experiences with workflows for automating data-intensive bioinformatics <https://biologydirect.biomedcentral.com/articles/10.1186/s13062-015-0071-8>`_. Biology Direct.
* Schramm et al. 2015. `Mutational dynamics between primary and relapse neuroblastomas <https://www.nature.com/ng/journal/v47/n8/full/ng.3349.html>`_. Nature Genetics.
* Berulava et al. 2015. `N6-Adenosine Methylation in MiRNAs <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118438>`_. PLOS ONE.
* The Genome of the Netherlands Consortium 2014. `Whole-genome sequence variation, population structure and demographic history of the Dutch population <https://www.nature.com/ng/journal/v46/n8/full/ng.3021.html>`_. Nature Genetics.
*  Patterson et al. 2014. `WhatsHap: Haplotype Assembly for Future-Generation Sequencing Reads <https://online.liebertpub.com/doi/10.1089/cmb.2014.0157>`_. Journal of Computational Biology.
* Fernández et al. 2014. `H3K4me1 marks DNA regions hypomethylated during aging in human stem and differentiated cells <https://genome.cshlp.org/content/25/1/27.long>`_. Genome Research.
* Köster et al. 2014. `Massively parallel read mapping on GPUs with the q-group index and PEANUT <https://peerj.com/articles/606/>`_. PeerJ.
* Chang et al. 2014. `TAIL-seq: Genome-wide Determination of Poly(A) Tail Length and 3′ End Modifications <https://www.cell.com/molecular-cell/abstract/S1097-2765(14)00121-X>`_. Molecular Cell.
* Althoff et al. 2013. `MiR-137 functions as a tumor suppressor in neuroblastoma by downregulating KDM1A <https://onlinelibrary.wiley.com/doi/10.1002/ijc.28091/abstract;jsessionid=33613A834E2A2FDCCA49246C23DF777E.f04t02>`_. International Journal of Cancer.
* Marschall et al. 2013. `MATE-CLEVER: Mendelian-Inheritance-Aware Discovery and Genotyping of Midsize and Long Indels <https://bioinformatics.oxfordjournals.org/content/29/24/3143.long>`_. Bioinformatics.
* Rahmann et al. 2013. `Identifying transcriptional miRNA biomarkers by integrating high-throughput sequencing and real-time PCR data <https://www.sciencedirect.com/science/article/pii/S1046202312002605>`_. Methods.
* Martin et al. 2013. `Exome sequencing identifies recurrent somatic mutations in EIF1AX and SF3B1 in uveal melanoma with disomy 3 <https://www.nature.com/ng/journal/v45/n8/full/ng.2674.html>`_. Nature Genetics.
* Czeschik et al. 2013. `Clinical and mutation data in 12 patients with the clinical diagnosis of Nager syndrome <https://link.springer.com/article/10.1007%2Fs00439-013-1295-2>`_. Human Genetics.
* Marschall et al. 2012. `CLEVER: Clique-Enumerating Variant Finder <https://bioinformatics.oxfordjournals.org/content/28/22/2875.long>`_. Bioinformatics.


.. toctree::
   :caption: Getting started
   :name: getting_started
   :hidden:
   :maxdepth: 1

   getting_started/installation
   tutorial/tutorial
   tutorial/short


.. toctree::
  :caption: Executing workflows
  :name: execution
  :hidden:
  :maxdepth: 1

  executing/cli
  executing/cluster-cloud
  executing/caching
  executing/interoperability

.. toctree::
    :caption: Defining workflows
    :name: snakefiles
    :hidden:
    :maxdepth: 1

    snakefiles/writing_snakefiles
    snakefiles/rules
    snakefiles/configuration
    snakefiles/modularization
    snakefiles/remote_files
    snakefiles/utils
    snakefiles/deployment
    snakefiles/reporting


.. toctree::
    :caption: API Reference
    :name: api-reference
    :hidden:
    :maxdepth: 1

    api_reference/snakemake
    api_reference/snakemake_utils
    api_reference/internal/modules


.. toctree::
    :caption: Project Info
    :name: project-info
    :hidden:
    :maxdepth: 1

    project_info/citations
    project_info/more_resources
    project_info/faq
    project_info/contributing
    project_info/authors
    project_info/history
    project_info/license
