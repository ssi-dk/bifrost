# Dockerfile to for:

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost__cge_mlst" \
    description="Docker environment for cge_mlst in bifrost" \
    version="1.0" \
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
    conda install -y biopython==1.74; \
    pip install bifrostlib==1.7; \
    pip install cgecore==1.5.0; \
    pip install tabulate==0.8.3; \
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
    git checkout f8ea197;
RUN \
    cd /bifrost_resources/resfinder_db; \
    python3 INSTALL.py kma_index; 