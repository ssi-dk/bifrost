# Dockerfile to for:

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost__ssi_stamper_check" \
    description="Docker environment for ssi_stamper in bifrost" \
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
    pip install -q bifrostlib==2.0.4; \
    wget -q https://raw.githubusercontent.com/ssi-dk/bifrost/master/setup/docker_umask_002.sh; \
    chmod +x  docker_umask_002.sh
ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]