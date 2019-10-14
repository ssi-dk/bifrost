FROM \
    ssidk/bifrost-base:2.0

LABEL \
    name="bifrost-min_read_check" \
    description="Docker environment for min_read_check in bifrost" \
    version="2.0" \
    DBversion="31/07/2019" \
    maintainer="kimn@ssi.dk;"

RUN \
    conda install -yq -c bioconda bbmap==38.58; \
    # created in base
    cd /bifrost_resources; \
    wget -q https://raw.githubusercontent.com/ssi-dk/bifrost/master/setup/adapters.fasta;

ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]