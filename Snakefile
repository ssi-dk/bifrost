#!/usr/bin/env python3
import re
import pandas
import sys
import os
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "scripts/"))
import serum  # all serum lib functions also have access to the config file to prevent reduce excess parameter passing
import pkg_resources

# config_file = pkg_resources.resource_filename(workflow.snakefile "/config/config.yaml")
configfile: os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")
# from pytools.persistent_dict import PersistentDict
# storage = PersistentDict("qcquickie_storage")
# snakemake -s ~/code/serumqc/snakefiles/serumqc.snake --config R1_reads={read_location} R2_reads={read_location} Sample=Test
# snakemake -s ~/code/serumqc/snakefiles/serumqc.snake --config R1_reads=~/test/data/nextseq/FHA3_S64_L555_R1_001.fastq.gz R2_reads=~/test/data/nextseq/FHA3_S64_L555_R2_001.fastq.gz Sample=Test
# requires --config R1_reads={read_location},R2_reads={read_location}
# snakemake -s ~/git.repositories/SerumQC-private/batch_run.snake --config run_folder=../../data/tiny/ sample_sheet=/srv/data/BIG/NGS_facility/assembly/2018/180117_NS500304_0140_N_WGS_91_AHWHHFAFXX/sample_sheet.xlsx

run_folder = str(config["run_folder"])
sample_sheet = str(config["sample_sheet"])
global_threads = config["global"]["threads"]
global_memory_in_GB = config["global"]["memory"]
# my understanding is all helps specify final output
onsuccess:
    print("Workflow complete")
    output = ["status.txt"]
    with open(output[0], "w") as status:
        status.write("Success")
onerror:
    print("Workflow error")
    output = ["status.txt"]
    with open(output[0], "w") as status:
        status.write("Failure")

rule all:
    input:
        "init_complete"

rule set_up_run:
    input:
        run_folder = run_folder,
        check_provided_sample_info = "check_provided_sample_info"
    output:
        samplesheet = "sample_sheet.xlsx",
        run_info_yaml = "init_complete"
    params:
        samplesheet = sample_sheet
    run:
        serum.check__run_folder(input.run_folder)
        serum.check__combine_sample_sheet_with_run_info(output.samplesheet)
        serum.initialize__run_from_run_info()
        serum.start_initialized_samples()
        serum.initialize_complete()
        # post steps


rule check__provided_sample_info:
    message:
        "Running step: {rule}"
    input:
        run_folder = run_folder
    output:
        samplesheet = "sample_sheet.xlsx",
        run_info_yaml = touch("check_provided_sample_info")
    params:
        samplesheet = sample_sheet
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        os.path.join(folder_name, "log/check__provided_sample_info.log")
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/check_provided_sample_info.py")


rule qcquickie_samples:
    input:
        run_config = "run.yaml",
        cmd_qcquickie = expand({sample}/cmd_qcquickie.sh, sample=config["samples"])