FROM continuumio/miniconda3:23.5.2-0

LABEL maintainer="Adam Hoffman <adamhoffman21@hotmail.ca>"
LABEL description="RNA-seq Differential Expression Analysis Pipeline"

# Set working directory
WORKDIR /workspace

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    curl \
    pigz \
    parallel \
    && rm -rf /var/lib/apt/lists/*

# Copy environment file
COPY environment.yml .

# Create conda environment
RUN conda env create -f environment.yml && \
    conda clean --all --yes

# Activate environment by default
ENV PATH /opt/conda/envs/rnaseq-pipeline/bin:$PATH
SHELL ["/bin/bash", "-c"]

# Verify installations
RUN conda run -n rnaseq-pipeline \
    snakemake --version && \
    STAR --version && \
    samtools --version && \
    featureCounts -v && \
    fastqc --version && \
    multiqc --version

# Copy pipeline files
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/
COPY notebooks/ notebooks/

# Create data directories
RUN mkdir -p data/raw data/reference results/{qc,alignment,counts,de_analysis,plots} logs

# Set entrypoint
ENTRYPOINT ["snakemake"]
CMD ["--help"]
