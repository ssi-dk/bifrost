# Dockerfile to build image

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost__assemblatron" \
    description="Docker environment for assemblatron in bifrost" \
    version="1.0" \
    DBversion="31/07/2019" \
    maintainer="kimn@ssi.dk;"

RUN \
    conda install conda==4.7.10; \
    conda config --add channels conda-forge; \
    conda config --add channels bioconda; \
    conda install prokka==1.13.7; \
    conda install pandas==0.25.0; \
    conda install bbmap==38.58; \
    conda install skesa==2.3.0; \
    conda install minimap2==2.17; \
    conda install samtools==1.9; \
    conda install cyvcf2==0.11.4; \
    pip install bifrostlib==1.4; \
    pip install quast==5.0.2; \
    conda install wget==1.20.1; \
    mkdir bifrost_resources;
RUN \
    cd /bifrost_resources; \
    wget https://raw.githubusercontent.com/ssi-dk/bifrost/dev/setup/adapters.fasta;