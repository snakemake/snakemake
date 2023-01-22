FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="d630028d41ac6564d03ea65c4edb2d271d33f36409fedab419872fe30d7df6b4"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: env.yaml
#   prefix: /conda-envs/_1219696ceec39be99ff366456467efb2
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - bcftools =1.10
RUN mkdir -p /conda-envs/_1219696ceec39be99ff366456467efb2
COPY env.yaml /conda-envs/_1219696ceec39be99ff366456467efb2/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/_1219696ceec39be99ff366456467efb2 --file /conda-envs/_1219696ceec39be99ff366456467efb2/environment.yaml && \
    mamba clean --all -y
