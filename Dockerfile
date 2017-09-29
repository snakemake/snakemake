FROM bitnami/minideb:jessie
MAINTAINER Johannes Köster <johannes.koester@tu-dortmund.de>
ADD . /tmp/repo
RUN install_packages wget
RUN wget -O- http://neuro.debian.net/lists/xenial.us-ca.full > /etc/apt/sources.list.d/neurodebian.sources.list
RUN wget -O- http://neuro.debian.net/_static/neuro.debian.net.asc | apt-key add -
RUN install_packages bzip2 ca-certificates singularity-container
RUN wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/conda/bin:${PATH}
ENV LANG C.UTF-8
ENV SHELL /bin/bash
RUN conda env update --name root --file /tmp/repo/environment.yml && conda clean --all -y
RUN pip install /tmp/repo
