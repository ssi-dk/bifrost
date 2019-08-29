# Dockerfile to for:

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost__min_read_check" \
    description="Docker environment for min_read_check in bifrost" \
    version="2.0" \
    DBversion="31/07/2019" \
    maintainer="kimn@ssi.dk;"

RUN \
    conda install conda==4.7.10; \
    conda config --add channels conda-forge; \
    conda config --add channels bioconda; \
    conda install -yq ruamel.yaml==0.16.5; \
    conda install -yq pymongo==3.8.0; \
    conda install -yq pandas==0.25.0; \
    conda install -yq python-magic==0.4.15; \
    conda install -yq bbmap==38.58; \
    conda install -yq wget==1.20.1; \
    pip install -q bifrostlib==2.0.1; \
    mkdir bifrost_resources; \
    cd /bifrost_resources; \
    wget https://raw.githubusercontent.com/ssi-dk/bifrost/dev/setup/adapters.fasta;
