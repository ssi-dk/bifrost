#!/usr/bin/env python3

# TODO: should refactor a lot of this Snakefile into a more snakemake orientated solution utilizing wildcards

import re
import sys
import os

configfile: os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")

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

if "group" in config:
    group = str(config["group"])
else:
    group = "NA"

partition = str(config["partition"])
global_threads = config["threads"]
global_memory_in_GB = config["memory"]
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
        group = group,
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
