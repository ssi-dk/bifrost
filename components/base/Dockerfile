FROM \
    continuumio/miniconda3:4.7.10

LABEL \
    name="bifrost-base" \
    description="Docker environment for common library and structure in bifrost" \
    version="2.0.7" \
    DBversion="14/10/2019" \
    maintainer="kimn@ssi.dk;"

RUN \
    pip install -q bifrostlib==2.0.7; \
    conda install -yq -c conda-forge -c bioconda -c default snakemake-minimal==5.7.1; \
    mkdir bifrost; \
    mkdir bifrost_resources;
