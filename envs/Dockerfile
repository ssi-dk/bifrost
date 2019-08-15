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
    conda install --name base --file /spec-file.txt; 
