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
    apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config; \
    export VERSION=1.11 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz ; \
    echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    . ~/.bashrc ; \
    export VERSION=3.3.0 && \
    mkdir -p $GOPATH/src/github.com/sylabs && \
    cd $GOPATH/src/github.com/sylabs && \
    wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz && \
    tar -xzf singularity-${VERSION}.tar.gz && \
    cd ./singularity && \
    ./mconfig && \
    make -C ./builddir && \
    make -C ./builddir install
