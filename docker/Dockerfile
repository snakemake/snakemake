# a docker image based on Ubuntu with snakemake installed
FROM ubuntu:14.04
MAINTAINER Johannes Köster <johannes.koester@tu-dortmund.de>
RUN apt-get -qq update
RUN apt-get install -qqy python3-setuptools python3-docutils python3-flask
RUN easy_install3 snakemake
ENTRYPOINT ["snakemake"]
