FROM mambaorg/micromamba
LABEL org.opencontainers.image.authors="Johannes Köster <johannes.koester@tu-dortmund.de>, Alp Tuna <alp.tuna@cern.ch>"

ENV LANG C.UTF-8
ENV SHELL /bin/bash
ENV APT_PKGS="bzip2 ca-certificates curl wget gnupg2 squashfs-tools git"

ARG UID=1000
ARG GID=1000

USER root

RUN apt-get update \
    && apt-get install -y --no-install-recommends ${APT_PKGS} \
    && apt-get clean \
    && rm -rf /var/lib/apt /var/lib/dpkg /var/lib/cache /var/lib/log \
    && micromamba create -q -y -c conda-forge -n apptainer apptainer \
    && micromamba create -q -y -c bioconda -c conda-forge -n snakemake \
        snakemake-minimal --only-deps \
    && eval "$(micromamba shell hook --shell bash)" \
    && micromamba activate /opt/conda/envs/snakemake \
    && micromamba install -c conda-forge conda \
    && micromamba clean --all -y \
    && groupadd -g $GID snakemake \
    && useradd -u $UID -g $GID -s /bin/bash snakemake \
    && chown -R snakemake:snakemake /opt/conda/envs

ENV PATH="/opt/conda/envs/snakemake/bin:/opt/conda/envs/apptainer/bin:${PATH}"
    
COPY . /tmp/repo
WORKDIR /tmp/repo
RUN chown -R snakemake:snakemake /tmp/repo

USER snakemake    

RUN pip install .[reports,messaging,pep]
