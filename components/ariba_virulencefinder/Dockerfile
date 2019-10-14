# Dockerfile to for:

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost__ariba_virulencefinder" \
    description="Docker environment for ariba_virulencefinder in bifrost" \
    version="2.0" \
    DBversion="21/08/19" \
    maintainer="kimn@ssi.dk;"

RUN \
    conda config --add channels conda-forge; \
    conda config --add channels bioconda; \
    conda install -yq ruamel.yaml==0.16.5 pymongo==3.8.0 pandas==0.25.0 python-magic==0.4.15 ariba==2.14.1; \
    pip install -q bifrostlib==2.0.4; \
    mkdir bifrost_resources; \
    cd /bifrost_resources; \
    wget -q https://raw.githubusercontent.com/ssi-dk/bifrost/master/setup/docker_umask_002.sh; \
    chmod +x  docker_umask_002.sh \
    mkdir virulencefinder; \
    cd /bifrost_resources/virulencefinder; \
    ariba getref virulencefinder virulencefinder --version 80c55fe; \
    ariba prepareref -f virulencefinder.fa -m virulencefinder.tsv ref_db;

ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]