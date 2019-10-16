FROM \
    ssidk/bifrost-base:2.0

LABEL \
    name="bifrost-ariba_plasmidfinder" \
    description="Docker environment for ariba_plasmidfinder in bifrost" \
    version="2.0" \
    DBversion="28/08/19" \
    maintainer="kimn@ssi.dk;"

RUN \
    # Next 3 are for make which is needed by ariba for install
    apt-get update -qq --fix-missing; \
    apt-get install -y -qq ariba=2.13.3+ds-1; \
    # In base image
    cd /bifrost_resources; \
    mkdir plasmidfinder; \
    cd /bifrost_resources/plasmidfinder; \
    ariba getref plasmidfinder plasmidfinder --version 1307168; \
    ariba prepareref -f plasmidfinder.fa -m plasmidfinder.tsv ref_db; 
ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]