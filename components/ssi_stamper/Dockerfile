# Dockerfile to for:

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost__ssi_stamper_check" \
    description="Docker environment for ssi_stamper in bifrost" \
    version="1.2" \
    DBversion="31/07/19" \
    maintainer="kimn@ssi.dk;"

RUN \
    conda install conda==4.7.10; \
    conda config --add channels conda-forge; \
    conda config --add channels bioconda; \
    pip install bifrostlib==1.3; 