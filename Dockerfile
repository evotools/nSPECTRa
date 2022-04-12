FROM continuumio/miniconda3 AS build

LABEL authors="andrea.talenti@ed.ac.uk" \
      description="Docker image containing base requirements for n-spectr pipelines"

# Install the package as normal:
COPY envs/all_environment.yml ./environment.yml

# Install mamba to speed up the process
RUN conda install -c conda-forge -y mamba

# Create the environment
RUN mamba env create -f environment.yml

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
FROM debian:buster AS runtime

# Install procps in debian to make it compatible with reporting
RUN apt-get update && apt install -y file procps g++ curl wget && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget https://github.com/marbl/meryl/releases/download/v1.3/meryl-1.3.Linux-amd64.tar.xz && \
  tar xJf meryl-1.3.Linux-amd64.tar.xz && \
  mv meryl-1.3/ meryl/ && \
  rm meryl-1.3.Linux-amd64.tar.xz

# Install datasets
ADD https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets /usr/local/bin/datasets 
RUN chmod a+x /usr/local/bin/datasets

# Copy /venv from the previous stage:
COPY --from=build /venv /venv

# When image is run, run the code with the environment
# activated:
ENV PATH /venv/bin/:/meryl/bin/:$PATH
SHELL ["/bin/bash", "-c"]
