FROM \
    continuumio/miniconda3:4.7.10

LABEL \
    name="bifrost-base" \
    description="Docker environment for common library and structure in bifrost" \
    version="2.1" \
    DBversion="14/10/2019" \
    maintainer="kimn@ssi.dk;"

RUN \
    pip install -q bifrostlib==2.0.4; \
    mkdir bifrost; \
    mkdir bifrost_resources; \
    cd bifrost_resources; \
    conda install -yq -c conda-forge -c bioconda -c default snakemake-minimal==5.7.1; \
    wget -q https://raw.githubusercontent.com/ssi-dk/bifrost/master/setup/docker_umask_002.sh; \
    chmod +x docker_umask_002.sh;

ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]