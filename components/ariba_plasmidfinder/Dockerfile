FROM \
    ssidk/bifrost-base:2.0

LABEL \
    name="bifrost-ariba_plasmidfinder" \
    description="Docker environment for ariba_plasmidfinder in bifrost" \
    version="2.0" \
    DBversion="28/08/19" \
    maintainer="kimn@ssi.dk;"

RUN \
    conda install -yq -c bioconda ariba==2.14.1; \
    # In base image
    cd /bifrost_resources; \
    mkdir plasmidfinder; \
    cd /bifrost_resources/plasmidfinder; \
    ariba getref plasmidfinder plasmidfinder --version 1307168; \
    ariba prepareref -f plasmidfinder.fa -m plasmidfinder.tsv ref_db; 

ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]