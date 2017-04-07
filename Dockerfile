FROM continuumio/miniconda3
MAINTAINER Johannes Köster <johannes.koester@tu-dortmund.de>
ADD . /tmp/repo
RUN apt-get upgrade -y
RUN conda update conda
RUN conda install conda-env-pip
RUN conda env update --name root --file /tmp/repo/environment.yml
RUN pip install /tmp/repo
