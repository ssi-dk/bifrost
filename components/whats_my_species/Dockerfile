# Dockerfile to for:

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost__whats_my_species_check" \
    description="Docker environment for whats_my_species in bifrost" \
    version="1.0" \
    DBversion="31/07/19" \
    maintainer="kimn@ssi.dk;"

RUN \
    conda install conda==4.7.10; \
    conda config --add channels conda-forge; \
    conda config --add channels bioconda; \
    conda install pandas==0.25.0; \
    conda install kraken==1.1.1; \
    conda install bracken==1.0.0; \
    pip install bifrostlib==1.2; \
    conda install wget==1.20.1; \
    # Directory for static DBs
    mkdir bifrost_resources;
RUN \
    # Kraken mini DB and bracken
    cd /bifrost_resources; \
    wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz; \
    # wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz; \
    tar -zxvf minikraken_20171019_4GB.tgz; \
    # tar -zxvf minikraken_20171019_8GB.tgz; \
    wget https://ccb.jhu.edu/software/bracken/dl/minikraken_4GB_100mers_distrib.txt;
    # wget https://ccb.jhu.edu/software/bracken/dl/minikraken_8GB_100mers_distrib.txt;

