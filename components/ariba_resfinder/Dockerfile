# Dockerfile to for:

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost__ariba_resfinder" \
    description="Docker environment for ariba_resfinder in bifrost" \
    version="2.0" \
    DBversion="21/08/19" \
    maintainer="kimn@ssi.dk;"

RUN \
    conda config --add channels conda-forge; \
    conda config --add channels bioconda; \
    conda install ariba==2.14.1; \
    pip install bifrostlib==2.0; \
    mkdir bifrost_resources; \
    cd /bifrost_resources; \
    mkdir resfinder; \
    cd /bifrost_resources/resfinder; \
    ariba getref resfinder resfinder --version f8ea197; \
