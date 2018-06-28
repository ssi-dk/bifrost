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
"serumqc"
# my understanding is all helps specify final output
onsuccess:
    print("Workflow complete")
    shell("touch serumqc_successfully_initialized_on_" + str(datetime.datetime.now()).replace(" ", "_"))

onerror:
    print("Workflow error")
    shell("touch serumqc_failed_to_initialized_on_" + str(datetime.datetime.now()).replace(" ", "_"))

ruleorder: setup > initialize_run


rule all:
    input:
        "run.yaml"


rule setup:
    output:
        directory = folder_name
    shell:
        "mkdir {output}"


rule get_git_hash_of_serumqc:
    input:
        "serumqc"
    output:
        "serumqc/git_hash.txt"
    shell:
        "git --git-dir {workflow.snakefile} rev-parse snakemake 1> {output}"

rule initialize_run:
    message:
        "Running step: {rule}"
    input:
        "serumqc/git_hash.txt",
        run_folder = run_folder,
    output:
        samplesheet = "sample_sheet.tsv",
        output = "run.yaml"
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

# can break this down to 2 parts where you create the sample_sheet in one and then prep for run with the other
# rule start_run:
#     input:
#         "run.yaml"
#     output:
#         touch("run_started")
#     shell:
#         "bash run_cmd_serumqc.sh"