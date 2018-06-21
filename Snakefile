#!/usr/bin/env python3

# TODO: should refactor a lot of this Snakefile into a more snakemake orientated solution utilizing wildcards

import re
import sys
import os
import datetime

configfile: os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")

components = str(config["components"])
run_folder = str(config["run_folder"])
sample_sheet = str(config["sample_sheet"])
group = str(config["group"])

partition = str(config["partition"])
global_threads = config["threads"]
global_memory_in_GB = config["memory"]
folder_name = "run_info"
# my understanding is all helps specify final output
onsuccess:
    print("Workflow complete")
    shell("touch qcquickie_successfully_initialized_on_" + str(datetime.datetime.now()).replace(" ", "_"))

onerror:
    print("Workflow error")
    shell("touch qcquickie_failed_to_initialized_on_" + str(datetime.datetime.now()).replace(" ", "_"))

ruleorder: setup > initialize_run


rule all:
    input:
        "run.yaml"


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
        config = config,
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

# rule check_provided_sample_sheet:

# rule 