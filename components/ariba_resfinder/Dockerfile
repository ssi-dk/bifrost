FROM \
    ssidk/bifrost-base:2.0

LABEL \
    name="bifrost-ariba_resfinder" \
    description="Docker environment for ariba_resfinder in bifrost" \
    version="2.0" \
    DBversion="21/08/19" \
    maintainer="kimn@ssi.dk;"

RUN \
    apt-get update -qq --fix-missing; \
    apt-get install -y -qq ariba=2.13.3+ds-1; \
    # In base image
    cd /bifrost_resources; \
    mkdir resfinder; \
    cd /bifrost_resources/resfinder; \
    ariba getref resfinder resfinder --version 149209d; \
    ariba prepareref -f resfinder.fa -m resfinder.tsv ref_db;

ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]