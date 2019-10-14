FROM \
    ssidk/bifrost-base:2.0

LABEL \
    name="bifrost-ariba_resfinder" \
    description="Docker environment for ariba_resfinder in bifrost" \
    version="2.0" \
    DBversion="21/08/19" \
    maintainer="kimn@ssi.dk;"

RUN \
    conda install -yq -c bioconda ariba==2.14.1; \
    # In base image
    cd /bifrost_resources; \
    mkdir resfinder; \
    cd /bifrost_resources/resfinder; \
    ariba getref resfinder resfinder --version 149209d; \
    ariba prepareref -f resfinder.fa -m resfinder.tsv ref_db;
