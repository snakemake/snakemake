ARG PYTHON_VERSION="3.8"
ARG CONDA_VERSION="latest"
ARG SINGULARITY_VERSION="3.6.1"
ARG NONROOT_USER="nonroot"
ARG NONROOT_UID="1000"
ARG NONROOT_GID="100"
ARG DEBIAN_FRONTEND="noninteractive"

# Compile Singularity.
FROM golang:1.14.6-stretch AS singularityBuilder
ARG DEBIAN_FRONTEND
ARG SINGULARITY_VERSION
WORKDIR /opt/singularity
RUN apt-get update -qq \
    && apt-get install --yes --no-install-recommends \
        cryptsetup \
        git \
        libgpgme11-dev \
        libseccomp-dev \
        libssl-dev \
        pkg-config \
        squashfs-tools \
        uuid-dev \
    && rm -rf /var/lib/apt/lists/* \
    && curl -fsSL "https://github.com/hpcng/singularity/releases/download/v${SINGULARITY_VERSION}/singularity-${SINGULARITY_VERSION}.tar.gz" \
    | tar xz --strip-components 1 \
    && ./mconfig --prefix=/opt/singularity \
    && cd builddir \
    && make \
    && make install

# Install Snakemake and its dependencies.
FROM "continuumio/miniconda3:${CONDA_VERSION}" as snakemakeBuilder
ARG DEBIAN_FRONTEND
ARG PYTHON_VERSION
RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    # Fast drop-in replacement for conda.
    && conda install --yes --quiet mamba \
    && mamba install --yes --quiet "python=$PYTHON_VERSION" \
    && mamba install --yes --quiet --only-deps snakemake snakemake-minimal \
    # Is this OK?
    && mamba remove --yes --quiet mamba \
    && conda clean --all --yes
WORKDIR /tmp/snakemake
COPY . .
RUN python -m pip install --no-cache-dir .[reports,messaging,google-cloud]

# Build final image.
FROM "continuumio/miniconda3:${CONDA_VERSION}"
ARG DEBIAN_FRONTEND
ARG NONROOT_USER
ARG NONROOT_UID
ARG NONROOT_GID
# Install Singularity runtime dependencies.
RUN apt-get update -qq \
    && apt-get install --yes --no-install-recommends \
        cryptsetup \
        squashfs-tools \
    && rm -rf /var/lib/apt/lists/*
# Create non-root user.
RUN useradd --create-home --shell /bin/bash -u "$NONROOT_UID" -g "$NONROOT_GID" "$NONROOT_USER"
USER "$NONROOT_UID"
# Add Singularity and conda. This overwrites the /opt/conda directory but that is OK
# because the snakemakeBuilder stage and the final stage are derived from the same
# base image.
COPY --from=snakemakeBuilder --chown="$NONROOT_UID":"$NONROOT_GID" /opt /opt
COPY --from=singularityBuilder /opt/singularity /opt/singularity
ENV PATH="/opt/conda/bin:/opt/singularity/bin:$PATH"
WORKDIR /work
ENTRYPOINT ["/opt/conda/bin/snakemake"]
LABEL maintainer="Johannes Köster <johannes.koester@tu-dortmund.de>"
