#!/usr/bin/env python3
import re
import sys
import os

# sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "scripts/"))
# import serum  # all serum lib functions also have access to the config file to prevent reduce excess parameter passing

# config_file = pkg_resources.resource_filename(workflow.snakefile "/config/config.yaml")
configfile: os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")
# from pytools.persistent_dict import PersistentDict
# storage = PersistentDict("qcquickie_storage")
# snakemake -s ~/code/serumqc/snakefiles/serumqc.snake --config R1_reads={read_location} R2_reads={read_location} Sample=Test
# snakemake -s ~/code/serumqc/snakefiles/serumqc.snake --config R1_reads=~/test/data/nextseq/FHA3_S64_L555_R1_001.fastq.gz R2_reads=~/test/data/nextseq/FHA3_S64_L555_R2_001.fastq.gz Sample=Test
# requires --config R1_reads={read_location},R2_reads={read_location}
# snakemake -s ~/git.repositories/SerumQC-private/batch_run.snake --config run_folder=../../data/tiny/ sample_sheet=/srv/data/BIG/NGS_facility/assembly/2018/180117_NS500304_0140_N_WGS_91_AHWHHFAFXX/sample_sheet.xlsx partition=daytime 

if "components" in config:
    components = str(config["components"])
else:
    components = "default"

if "run_folder" in config:
    run_folder = str(config["run_folder"])
elif "samples" in os.getcwd():
    run_folder = "samples"

if "sample_sheet" in config:
    sample_sheet = str(config["sample_sheet"])
else:
    sample_sheet = ""

partition = str(config["partition"])
global_threads = config["global"]["threads"]
global_memory_in_GB = config["global"]["memory"]
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

ruleorder: setup > initialize_run > print_run

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
        run_folder = run_folder,
    output:
        samplesheet = "sample_sheet.tsv",
        run_config_yaml = "run.yaml",
        check = touch(os.path.join(folder_name, "initialize_run"))
    params:
        samplesheet = sample_sheet,
        partition = partition,
        components = components,
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "envs/python_packages.yaml"
    log:
        os.path.join(folder_name, "log/initialize_run.log")
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/initialize_run.py")


def get_sample_names(wildcards):
    return [directory for directory in os.listdir() if os.path.isdir(directory) and os.path.isfile(os.path.join(directory, "sample.yaml"))]


rule print_run:
    message:
        "Running step: {rule}"
    input:
        samples = get_sample_names,
        run_config = "run.yaml"
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
        print("Found samples: {}".format(input.samples))
