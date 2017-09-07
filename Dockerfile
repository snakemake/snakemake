FROM bitnami/minideb:jessie
MAINTAINER Johannes Köster <johannes.koester@tu-dortmund.de>
ADD . /tmp/repo
# taken from condaforge/linux-anvil 
#RUN apt-get update && \
#    apt-get install -y wget bzip2 && \
#    rm -rf /var/lib/apt/lists/*
RUN install_packages wget bzip2
RUN wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/conda/bin:${PATH}
ENV LANG C.UTF-8
ENV SHELL /bin/bash
RUN conda env update --name root --file /tmp/repo/environment.yml && conda clean --all -y
RUN pip install /tmp/repo
