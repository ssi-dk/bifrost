#!/usr/bin/env python3
import re
import pandas
import sys
import os
from ruamel.yaml import YAML
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "scripts/"))
import serum  # all serum lib functions also have access to the config file to prevent reduce excess parameter passing

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
global run_config
folder_name = "run_info"
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

ruleorder: setup > initialize_run > set_run_info > print_run

rule all:
    input:
        os.path.join(folder_name, "print_run")


rule setup:
    output:
        dir = folder_name
    shell:
        "mkdir {output}"


rule initialize_run:
    message:
        "Running step: {rule}"
    input:
        run_folder = run_folder
    output:
        samplesheet = "sample_sheet.xlsx",
        run_config = "run.yaml",
        check = touch(os.path.join(folder_name, "initialize_run"))
    params:
        samplesheet = sample_sheet
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        os.path.join(folder_name, "log/initialize_run.log")
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/initialize_run.py")


rule set_run_info:
    message:
        "Running step: {rule}"
    input:
        run_config = "run.yaml"
    output:
        check = touch(os.path.join(folder_name, "set_run_info"))
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        os.path.join(folder_name, "log/set_run_info.log")
    run:
        yaml = YAML(typ='safe')
        yaml.default_flow_style = False
        with open(input.run_config, "r") as yaml_stream:
            config = yaml.load(yaml_stream)

rule print_run:
    message:
        "Running step: {rule}"
    input:
        os.path.join(folder_name, "set_run_info")
    output:
        check = touch(os.path.join(folder_name, "print_run"))
    params:
        samplesheet = sample_sheet
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        os.path.join(folder_name, "log/print_run.log")
    run:
        print(run_config)
