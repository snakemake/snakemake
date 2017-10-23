FROM bitnami/minideb:stretch
MAINTAINER Johannes Köster <johannes.koester@tu-dortmund.de>
ENV SINGULARITY_VERSION=2.3.2
ADD . /tmp/repo
RUN install_packages wget bzip2 ca-certificates gnupg2
RUN wget -O- http://neuro.debian.net/lists/xenial.us-ca.full > /etc/apt/sources.list.d/neurodebian.sources.list
RUN wget -O- http://neuro.debian.net/_static/neuro.debian.net.asc | apt-key add -
RUN install_packages singularity-container
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/conda/bin:${PATH}
ENV LANG C.UTF-8
ENV SHELL /bin/bash
RUN conda env update --name root --file /tmp/repo/environment.yml && conda clean --all -y
RUN pip install /tmp/repo
