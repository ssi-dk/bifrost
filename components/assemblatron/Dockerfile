FROM \
    ssidk/bifrost-base:2.0

LABEL \
    name="bifrost-assemblatron" \
    description="Docker environment for assemblatron in bifrost" \
    version="2.0" \
    DBversion="31/07/2019" \
    maintainer="kimn@ssi.dk;"

RUN \
    conda install -yq -c bioconda prokka==1.13.7; \
    conda install -yq -c bioconda bbmap==38.58; \
    conda install -yq -c bioconda skesa==2.3.0; \
    conda install -yq -c bioconda minimap2==2.17; \
    conda install -yq -c bioconda samtools==1.9; \
    conda install -yq -c bioconda cyvcf2==0.11.4; \
    pip install -q quast==5.0.2; \
    # In base image
    cd /bifrost_resources; \
    wget -q https://raw.githubusercontent.com/ssi-dk/bifrost/master/setup/adapters.fasta;

ENTRYPOINT \
    ["/bifrost_resources/docker_umask_002.sh"]