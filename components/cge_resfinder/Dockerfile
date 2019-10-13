# Dockerfile to for:

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost__cge_resfinder" \
    description="Docker environment for cge_resfinder in bifrost" \
    version="2.0" \
    DBversion="19/07/2019" \
    maintainer="kimn@ssi.dk;"

RUN \
    # For 'make' needed for kma
    apt-get update -qq --fix-missing; \
    apt-get install -y -qq build-essential; \
    apt-get install -y -qq zlib1g-dev; \
    conda install -y conda==4.7.10; \
    conda config --add channels conda-forge; \
    conda config --add channels bioconda; \
    conda install -yq ruamel.yaml==0.16.5; \
    conda install -yq biopython==1.74.0; \
    conda install -yq pymongo==3.8.0; \
    conda install -yq pandas==0.25.0; \
    conda install -yq python-magic==0.4.15; \
    pip install -q bifrostlib==2.0.4; \
    pip install -q cgecore==1.5.0; \
    pip install -q tabulate==0.8.3; \
    # Directory for static DBs
    mkdir bifrost_resources; \
    cd /bifrost_resources; \
    # Install KMA
    git clone --branch 1.2.10a https://bitbucket.org/genomicepidemiology/kma.git; \
    cd kma; \
    make;
ENV PATH $PATH":/bifrost_resources/kma/"
RUN \
    cd /bifrost_resources; \
    # Install cge Resfinder
    git clone https://bitbucket.org/genomicepidemiology/resfinder.git; \
    cd resfinder; \
    git checkout d210e15;
ENV PATH $PATH":/bifrost_resources/resfinder/"
RUN \
    cd /bifrost_resources; \
    # Install cge Resfinder DB
    # Grabbing master as the versions aren't tagged
    git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git; \
    # grab specific commit for date assurance 2019/07/19
    cd resfinder_db; \ 
    git checkout 149209d;
RUN \
    cd /bifrost_resources/resfinder_db; \
    python3 INSTALL.py kma_index; \
    wget -q https://raw.githubusercontent.com/ssi-dk/bifrost/master/setup/docker_umask_002.sh; \
    chmod +x  docker_umask_002.sh
ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]