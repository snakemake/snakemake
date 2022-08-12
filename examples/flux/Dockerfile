FROM fluxrm/flux-sched:focal
# docker build -t flux-snake .
# docker run -it flux-snake
USER root
#ENV PATH /home/fluxuser/.local/bin:${PATH}
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=/opt/conda:/bin:$PATH
RUN pip install git+https://github.com/vsoch/snakemake@add/flux-executor
USER fluxuser
COPY ./Snakefile /home/fluxuser/Snakefile
