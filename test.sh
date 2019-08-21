#!/bin/bash
echo $BIFROST_DB_PATH > /bifrost_key.txt;
echo $BIFROST_DB_PATH

export BIFROST_DB_KEY=/bifrost_key.txt;


mkdir -p /bifrost_run/samples; 
cd /bifrost_run/samples; 
wget -O SRR801237_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR801/SRR801237/SRR801237_1.fastq.gz; 
wget -O SRR801237_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR801/SRR801237/SRR801237_2.fastq.gz; 
wget -O SRR801202_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR801/SRR801202/SRR801202_1.fastq.gz; 
wget -O SRR801202_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR801/SRR801202/SRR801202_2.fastq.gz; 

cd /bifrost_run; 
snakemake -s src/bifrost.smk --config \
    read_pattern="(?P<sample_name>.+?)_(?P<paired_read_number>R[1|2])(?P<file_extension>\.fastq\.gz)" \
    run_name=test_run \
    grid=none \
    singularity_prefix="/singularity" \
    components=min_read_check;

cat /bifrost_run/bifrost/log/*;