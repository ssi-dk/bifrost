# Dockerfile to build image

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost__assemblatron" \
    description="Docker environment for assemblatron in bifrost" \
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
    conda install -yq prokka==1.13.7; \
    conda install -yq bbmap==38.58; \
    conda install -yq skesa==2.3.0; \
    conda install -yq minimap2==2.17; \
    conda install -yq samtools==1.9; \
    conda install -yq cyvcf2==0.11.4; \
    pip install -q bifrostlib==2.0.4; \
    pip install -q quast==5.0.2; \
    conda install wget==1.20.1; \
    mkdir bifrost_resources; \
    cd /bifrost_resources; \
    wget -q https://raw.githubusercontent.com/ssi-dk/bifrost/master/setup/docker_umask_002.sh; \
    chmod +x  docker_umask_002.sh \
    wget -q https://raw.githubusercontent.com/ssi-dk/bifrost/master/setup/adapters.fasta;

ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]