FROM \
    ssidk/bifrost-base:2.0

LABEL \
    name="bifrost-cge_resfinder" \
    description="Docker environment for cge_resfinder in bifrost" \
    version="2.0" \
    DBversion="19/07/2019" \
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
    python3 INSTALL.py kma_index;

ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]