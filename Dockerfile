FROM continuumio/anaconda
MAINTAINER Johannes Köster <johannes.koester@tu-dortmund.de>
ADD environment.yml /tmp/environment.yml
RUN conda update conda-env
RUN conda env update --name root --file /tmp/environment.yml
RUN conda install -c bioconda snakemake
ENTRYPOINT ["/bin/bash"]
