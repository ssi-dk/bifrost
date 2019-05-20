#!/bin/sh
# Assuming installed mongodb and some set up for it

conda env create --file envs/bifrost_for_install_full.yaml
source activate bifrost;

mkdir resources;
cd resources;
wget http://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz;
tar -xzf minikraken_20171019_8GB.tgz;
cd minikraken_20171019_8GB
wget https://ccb.jhu.edu/software/bracken/dl/minikraken_8GB_100mers_distrib.txt;
cd ..
ariba getref resfinder ariba_resfinder_db


mkdir tmp_programs;
cd tmp_programs;
git clone git@bitbucket.org:genomicepidemiology/kma.git .;
git clone git@bitbucket.org:genomicepidemiology/mlst.git .;
