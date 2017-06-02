FROM debian
MAINTAINER Johannes Köster <johannes.koester@tu-dortmund.de>
ADD . /tmp/repo
# taken from condaforge/linux-anvil 
RUN curl -s -L https://repo.continuum.io/miniconda/Miniconda3-4.2.12-Linux-x86_64.sh > miniconda.sh && \
    openssl md5 miniconda.sh | grep d0c7c71cc5659e54ab51f2005a8d96f3 && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    export PATH=/opt/conda/bin:$PATH && \
    conda update --all --yes && \
    conda clean -tipy
RUN conda env update --name root --file /tmp/repo/environment.yml
RUN pip install /tmp/repo
