FROM condaforge/miniforge3:26.1.1-3 AS build

LABEL authors="andrea.talenti@ed.ac.uk" \
      description="Docker image containing base requirements for n-spectra pipelines"

# Install compiler
RUN apt update -qq && apt install -y -qq build-essential

# Create the environment
RUN mamba create -n nspectra -c conda-forge -c bioconda -y python=3.10 \
    pip matplotlib numpy>=1.22.4 pandas polars==1.9.0 \
    pyarrow pysam scipy bioconda::java-jdk
RUN mamba install -n nspectra -y \
    bioconda::samtools bioconda::bcftools bioconda::vcftools \
    bioconda::bedtools bioconda::plink=1.90
RUN mamba install -n nspectra -y bioconda::ucsc-twobitinfo bioconda::ucsc-fatotwobit \
    bioconda::ucsc-wigtobigwig bioconda::ucsc-bigwigtobedgraph
RUN mamba install -n nspectra -y bioconda::shapeit4
RUN mamba install -n nspectra -y bioconda::shapeit5
RUN mamba install -n nspectra -y \
    bioconda::perl-bioperl \
    bioconda::phast=1.9.9 bioconda::mutyper bioconda::perl-bio-db-hts \
    bioconda::tabix bioconda::vcflib
RUN mamba install -n nspectra -y \
    conda-forge::r-base>=4.1.0 \
    conda-forge::r-cowplot \
    conda-forge::r-tidyverse \
    conda-forge::r-reshape2 \
    conda-forge::r-ggfortify \
    conda-forge::r-ggforce \
    conda-forge::r-lazyeval \
    conda-forge::r-plyr

# Install conda-pack:
RUN mamba install -c conda-forge conda-pack

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack -n nspectra -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /venv/bin/conda-unpack


# The runtime-stage image; we can use Debian as the
# base image since the Conda env also includes Python
# for us.
FROM ubuntu:24.04 AS runtime

# Install procps in debian to make it compatible with reporting
RUN apt-get update && apt install -y file procps g++ curl git wget parallel && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install datasets
ADD https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets /usr/local/bin/datasets 
RUN chmod a+x /usr/local/bin/datasets

# Copy /venv from the previous stage:
COPY --from=build /venv /venv

# When image is run, run the code with the environment
# activated:
ENV PATH=/venv/bin/:$PATH
SHELL ["/bin/bash", "-c"]
