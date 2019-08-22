# Dockerfile to for:

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost__cge_mlst" \
    description="Docker environment for cge_mlst in bifrost" \
    version="2.0" \
    DBversion="29/07/2019" \
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
    pip install bifrostlib==2.0; \
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
    # Install cge MLST
    # Grabbing master as I need a version which is untagged
    git clone https://bitbucket.org/genomicepidemiology/mlst.git; \
    cd mlst; \
    # Grabbing specific version with matrix output
    git checkout 96e62e3;
ENV PATH $PATH":/bifrost_resources/mlst/"
RUN \
    cd /bifrost_resources; \
    # Install cge MLST DB
    # Grabbing master as the versions aren't tagged
    git clone https://bitbucket.org/genomicepidemiology/mlst_db.git; \
    # grab specific commit for date assurance 2019/07/29
    cd mlst_db; \
    git checkout ced7b97;
ENV MLST_DB "/bifrost_resources/mlst_db" 
RUN \
    cd /bifrost_resources/mlst_db; \
    python3 INSTALL.py kma_index; 