# Dockerfile to create the bifrost environment image

# Load miniconda Docker image to work off of
FROM continuumio/miniconda3:4.6.14

LABEL \
    name="bifrost-env" \
    description="Docker environment for bifrost" \
    version="1.0" \
    maintainer="mbas@ssi.dk;"

RUN mkdir /bifrost_run

# Copying bifrost info
COPY spec-file.txt /spec-file.txt

# Install requirements for bifrost
RUN \
    apt-get update && apt-get install -y squashfs-tools ; \
    conda install --name base --file /spec-file.txt; \
    ln -s /opt/conda/bin/python /usr/bin/python ; \
    wget -O- http://neuro.debian.net/lists/stretch.de-m.full | tee /etc/apt/sources.list.d/neurodebian.sources.list; \
    apt-get install -y gnupg2; \
    wget http://neuro.debian.net/_static/neuro.debian.net.asc ; \
    apt-key add neuro.debian.net.asc ; \
    apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9; \
    apt-get update; \
    apt-get install -y singularity-container
