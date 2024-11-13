FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="cb075ce12360612250e0065299b188ea3547727b9abc02b9927671af50d26bca"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: env.yaml
#   prefix: /conda-envs/2a132eb0d84044eaf02032bd9f9851fd
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bcftools =1.21
RUN mkdir -p /conda-envs/2a132eb0d84044eaf02032bd9f9851fd
COPY env.yaml /conda-envs/2a132eb0d84044eaf02032bd9f9851fd/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/2a132eb0d84044eaf02032bd9f9851fd --file /conda-envs/2a132eb0d84044eaf02032bd9f9851fd/environment.yaml && \
    conda clean --all -y
