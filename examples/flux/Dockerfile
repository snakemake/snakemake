FROM fluxrm/flux-sched:focal
# docker build -t flux-snake .
# docker run -it flux-snake
USER root
ENV PATH=/opt/conda:/bin:$PATH
RUN pip install git+https://github.com/snakemake/snakemake@main && \
    ln -s /bin/python3 /usr/local/bin/python
USER fluxuser
COPY ./Snakefile /home/fluxuser/Snakefile
