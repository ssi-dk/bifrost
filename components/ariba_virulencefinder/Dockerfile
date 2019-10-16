FROM \
    ssidk/bifrost-base:2.0

LABEL \
    name="bifrost-ariba_virulencefinder" \
    description="Docker environment for ariba_virulencefinder in bifrost" \
    version="2.0" \
    DBversion="21/08/19" \
    maintainer="kimn@ssi.dk;"

RUN \
    apt-get update -qq --fix-missing; \
    apt-get install -y -qq ariba=2.13.3+ds-1; \
    # In base image
    cd /bifrost_resources; \
    mkdir virulencefinder; \
    cd /bifrost_resources/virulencefinder; \
    ariba getref virulencefinder virulencefinder --version 80c55fe; \
    ariba prepareref -f virulencefinder.fa -m virulencefinder.tsv ref_db;

ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]