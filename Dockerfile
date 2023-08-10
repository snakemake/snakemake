FROM mambaorg/micromamba
LABEL org.opencontainers.image.authors="Johannes Köster <johannes.koester@tu-dortmund.de>"
ADD . /tmp/repo
WORKDIR /tmp/repo
ENV LANG C.UTF-8
ENV SHELL /bin/bash
USER root 

RUN micromamba create -q -y -c bioconda -c conda-forge -n snakemake \
    snakemake snakemake-minimal --only-deps && \
    eval "$(micromamba shell hook --shell bash)" && \
    micromamba activate /opt/conda/envs/snakemake && \
    micromamba install -c conda-forge singularity && \
    micromamba clean --all -y 

ENV PATH /opt/conda/envs/snakemake/bin:${PATH}
RUN pip install .[reports,messaging,google-cloud,azure]