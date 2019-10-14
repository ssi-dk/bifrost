# Dockerfile to for:

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost__whats_my_species_check" \
    description="Docker environment for whats_my_species in bifrost" \
    version="2.0" \
    DBversion="31/07/19" \
    maintainer="kimn@ssi.dk;"

RUN \
    conda install conda==4.7.10; \
    conda config --add channels conda-forge; \
    conda config --add channels bioconda; \
    conda install -yq ruamel.yaml==0.16.5; \
    conda install -yq pymongo==3.8.0; \
    conda install -yq pandas==0.25.0; \
    conda install -yq python-magic==0.4.15; \
    conda install -yq kraken==1.1.1; \
    conda install -yq bracken==1.0.0; \
    conda install -yq wget==1.20.1; \
    pip install -q bifrostlib==2.0.4; \
    # Directory for static DBs
    mkdir bifrost_resources;
RUN \
    # Kraken mini DB and bracken
    cd /bifrost_resources; \
    wget -q https://raw.githubusercontent.com/ssi-dk/bifrost/master/setup/docker_umask_002.sh; \
    chmod +x  docker_umask_002.sh \
    # NOTE: even though the file is dated 20171019 it's actually when unzipped 20171013, to avoid some issues im renaming it minikraken
    # wget -q https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz; \
    wget -q https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz; \
    mkdir minikraken; \
    # tar -zxf minikraken_20171019_4GB.tgz -C minikraken --strip-components=1; \
    tar -zxvf minikraken_20171019_8GB.tgz; \
    rm minikraken_20171019_4GB.tgz; \
    # wget -q https://ccb.jhu.edu/software/bracken/dl/minikraken_4GB_100mers_distrib.txt;
    wget -q https://ccb.jhu.edu/software/bracken/dl/minikraken_8GB_100mers_distrib.txt;

ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]
