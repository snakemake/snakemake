FROM bitnami/minideb:buster
MAINTAINER Johannes Köster <johannes.koester@tu-dortmund.de>
ADD . /tmp/repo
WORKDIR /tmp/repo
ENV PATH /opt/conda/bin:${PATH}
ENV LANG C.UTF-8
ENV SHELL /bin/bash
RUN install_packages wget curl bzip2 ca-certificates gnupg2 squashfs-tools git
RUN /bin/bash -c "curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh > mambaforge.sh && \
    bash mambaforge.sh -b -p /opt/conda && \
    conda config --system --set channel_priority strict && \
    rm mambaforge.sh"
RUN /bin/bash -c "mamba create -q -y -c conda-forge -c bioconda -n snakemake snakemake snakemake-minimal --only-deps && \
    source activate snakemake && \
    mamba install -q -y -c conda-forge singularity && \
    conda clean --all -y && \
    which python && \
    pip install .[reports,messaging,google-cloud]"
RUN echo "source activate snakemake" > ~/.bashrc
ENV PATH /opt/conda/envs/snakemake/bin:${PATH}
