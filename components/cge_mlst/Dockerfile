FROM \
    ssidk/bifrost-base:2.0

LABEL \
    name="bifrost-cge_mlst" \
    description="Docker environment for cge_mlst in bifrost" \
    version="2.0" \
    DBversion="29/07/2019" \
    maintainer="kimn@ssi.dk;"

RUN \
    # For 'make' needed for kma
    apt-get update -qq --fix-missing; \
    apt-get install -y -qq build-essential; \
    apt-get install -y -qq zlib1g-dev; \
    pip install -q cgecore==1.5.0; \
    pip install -q tabulate==0.8.3; \
    pip install -q biopython==1.74; \
    # In base image
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

ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]