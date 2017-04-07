FROM continuumio/miniconda3
MAINTAINER Johannes Köster <johannes.koester@tu-dortmund.de>
ADD . /tmp/repo
RUN apt-get upgrade -y
RUN conda update --quiet conda
RUN conda env update --name root --quiet --file /tmp/repo/environment.yml
RUN pip install /tmp/repo
