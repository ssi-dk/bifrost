#!/bin/sh

conda install -y create --file ../envs/bifrost_for_install_full.yaml; \
#Currently a bug with samtools installation, need to redo it from this source
conda install -y samtools conda-forge::ncurses; \
conda activate bifrost; \
#pip install bifrostlib; \
cd ../../; \
mkdir bifrost_resources; \
mv install/adapter.fasta bifrost_resources; \
cd bifrost_resources; \
wget http://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz; \
tar -xzf minikraken_20171019_8GB.tgz; \
rm minikraken_20171019_8GB.tgz; \
cd minikraken_20171019_8GB; \
wget https://ccb.jhu.edu/software/bracken/dl/minikraken_8GB_100mers_distrib.txt; \
cd ..; \
cd bifrost; \
ln -s ../bifrost_resources resources; \
