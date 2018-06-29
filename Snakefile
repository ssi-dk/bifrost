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
        "serumqc_setup_complete"


rule setup:
    output:
        directory = "serumqc"
    shell:
        "mkdir {output}"


rule initialize_run:
    message:
        "Running step: {rule}"
    input:
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
        "serumqc/log/initialize_run.log"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/initialize_run.py")


rule get_git_hash_of_serumqc:
    input:
        "run.yaml"
    output:
        "serumqc/git_hash.txt"
    run:
        with open(input.run_info_yaml_path, "r") as run_info_yaml:
            run_info = yaml.load(run_info_yaml)
        shell("git --git-dir {workflow.basedir}/.git rev-parse snakemake 1> {output}")
        with open(output, "r") as git_info:
            git_hash = git_info.readlines()[0].strip()
        run_info["run"]["git_hash"] = git_hash
        with open(input.read_info_yaml_path, "w") as run_info_yaml:
            yaml.dump(run_info, run_info_yaml)


rule get_conda_env:
    input:
        git_hash = "serumqc/git_hash.txt",
        run_info_yaml_path = "run.yaml"
    output:
        "serumqc/conda.yaml"
    run:
        with open(input.run_info_yaml_path, "r") as run_info_yaml:
            run_info = yaml.load(run_info_yaml)
        shell("conda env export 1> {output}")
        with open(output, "r") as conda_info_yaml:
            conda_info = yaml.load(conda_info_yaml)
        run_info["run"]["conda_env"] = conda_info
        with open(input.read_info_yaml_path, "w") as run_info_yaml:
            yaml.dump(run_info, run_info_yaml)

rule create_end_file:
    input:
        "serumqc/conda.yaml"
    output:
        "serumqc_setup_complete"
    shell:
        "touch {output}"

# can break this down to 2 parts where you create the sample_sheet in one and then prep for run with the other
# rule start_run:
#     input:
#         "run.yaml"
#     output:
#         touch("run_started")
#     shell:
#         "bash run_cmd_serumqc.sh"