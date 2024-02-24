.. _manual-main:

=========
Snakemake
=========

.. image:: https://img.shields.io/badge/Gitpod-ready--to--code-blue?color=%23022c22
    :target: https://gitpod.io/#https://github.com/snakemake/snakemake

.. image:: https://img.shields.io/conda/dn/bioconda/snakemake.svg?label=Bioconda&color=%23064e3b
    :target: https://bioconda.github.io/recipes/snakemake/README.html

.. image:: https://img.shields.io/pypi/pyversions/snakemake.svg?color=%23065f46
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/snakemake.svg?color=%23047857
    :target: https://pypi.python.org/pypi/snakemake

.. image:: https://img.shields.io/github/actions/workflow/status/snakemake/snakemake/docker-publish.yml?label=docker%20container&branch=main&color=%23059669
    :target: https://hub.docker.com/r/snakemake/snakemake

.. image:: https://img.shields.io/github/actions/workflow/status/snakemake/snakemake/main.yml?label=tests&color=%2310b981
    :target: https://github.com/snakemake/snakemake/actions?query=branch%3Amain+workflow%3ACI

.. image:: https://img.shields.io/badge/stack-overflow-orange.svg?color=%2334d399
    :target: https://stackoverflow.com/questions/tagged/snakemake

.. image:: https://img.shields.io/discord/753690260830945390?label=discord%20chat&color=%23a7f3d0
    :alt: Discord
    :target: https://discord.gg/NUdMtmr

.. image:: https://img.shields.io/twitter/follow/johanneskoester.svg?style=social&label=Follow&color=%236ee7b7
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

Snakemake is **highly popular**, with `>7 new citations per week <https://badge.dimensions.ai/details/id/pub.1018944052>`_.
For an introduction, please visit https://snakemake.github.io.

This is the home of the **Snakemake API documentation**, which is meant **purely for developers of snakemake's internal functionalities**.
If you are looking for the Snakemake user documentation (for writing workflows or executing them), please visit https://snakemake.readthedocs.io.



.. toctree::
    :caption: API Reference
    :name: api-reference
    :hidden:
    :maxdepth: 1

    api_reference/snakemake_api
    api_reference/snakemake_utils
    api_reference/internal/modules