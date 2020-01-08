FROM ssidk/bifrost-base:2.0.5

ARG version="2.0.5"
ARG last_updated="19/07/2019"
ARG name="cge_mlst"
ARG full_name="bifrost-${name}"

LABEL \
    name=${name} \
    description="Docker environment for ${full_name}" \
    version=${version} \
    resource_version=${last_updated} \
    maintainer="kimn@ssi.dk;"

#- Tools to install:start---------------------------------------------------------------------------
RUN \
    # For 'make' needed for kma
    apt-get update && apt-get install -y -qq --fix-missing \
        build-essential \
        zlib1g-dev; \
    pip install -q \
        cgecore==1.5.0 \
        tabulate==0.8.3 \
        biopython==1.74;
# KMA
RUN cd /bifrost_resources && \
    git clone --branch 1.2.10a https://bitbucket.org/genomicepidemiology/kma.git && \
    cd kma && \
    make;
ENV PATH /bifrost_resources/kma/:$PATH
# Resfinder
RUN cd /bifrost_resources && \
    git clone https://bitbucket.org/genomicepidemiology/mlst.git && \
    cd mlst && \
    git checkout 96e62e3;
ENV PATH /bifrost_resources/mlst/:$PATH
#- Tools to install:end ----------------------------------------------------------------------------

#- Additional resources (files/DBs): start ---------------------------------------------------------
# Resfinder DB from 2019/07/19
RUN cd /bifrost_resources && \
    git clone https://git@bitbucket.org/genomicepidemiology/mlst_db.git && \
    cd mlst_db && \ 
    git checkout ced7b97 && \
    python3 INSTALL.py kma_index;
#- Additional resources (files/DBs): end -----------------------------------------------------------

#- Source code:start -------------------------------------------------------------------------------
RUN cd /bifrost && \
    git clone --branch ${version} https://github.com/ssi-dk/${full_name}.git ${name};
#- Source code:end ---------------------------------------------------------------------------------

#- Set up entry point:start ------------------------------------------------------------------------
ENV PATH /bifrost/${name}/:$PATH
ENTRYPOINT ["launcher.py"]
CMD ["launcher.py", "--help"]
#- Set up entry point:end --------------------------------------------------------------------------