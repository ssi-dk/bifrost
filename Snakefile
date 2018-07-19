#!/usr/bin/env python3

# TODO: should refactor a lot of this Snakefile into a more snakemake orientated solution utilizing wildcards
import re
import sys
import os
import datetime
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "scripts"))
import datahandling

configfile: os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")

#Saving the config
component = "serumqc"

datahandling.save_yaml(config, "serumqc_config.yaml")

components = config["components"]
run_folder = config["run_folder"]
sample_sheet = config["sample_sheet"]
group = config["group"]
partition = config["partition"]
global_threads = config["threads"]
global_memory_in_GB = config["memory"]


onsuccess:
    print("Workflow complete")
    shell("touch serumqc_successfully_initialized_on_" + str(datetime.datetime.now()).replace(" ", "_"))


onerror:
    print("Workflow error")
    shell("touch serumqc_failed_to_initialized_on_" + str(datetime.datetime.now()).replace(" ", "_"))


rule all:
    input:
        component + "/" + component + "_complete"


rule setup:
    output:
        folder = directory(component)
    shell:
        "mkdir {output}"


rule_name = "generate_git_hash"
rule generate_git_hash:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        component
    output:
        component + "/git_hash.txt"
    shell:
        "git --git-dir {workflow.basedir}/.git rev-parse snakemake 1> {output} 2> {log.err_file}"


rule_name = "export_conda_env"
rule export_conda_env:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        component
    output:
        component + "/conda_env.yaml"
    shell:
        "conda env export 1> {output} 2> {log.err_file}"


rule initialize_components:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        git_hash = component + "/git_hash.txt",
        conda_env = component + "/conda_env.yaml"
    output:
        touch(component + "/initialize_components_complete"),
    run:
        sys.stdout.write("Started initialize_components\n")
        component_info = {}
        with open(input.git_hash, "r") as git_info:
            git_hash = git_info.readlines()[0].strip()
            component_info["git_hash"] = git_hash

        component_info["conda_env"] = datahandling.load_yaml(input.conda_env)
        component_info["config"]: config

        for component_name in components:
            component_info["name"] = component_name
            datahandling.save_component(component_info, component + "/" + component_name + ".yaml")
        shell("touch initialize_components_complete")
        print(component_info)
        sys.stdout.write("Done initialize_components\n")

# rule_name = "initialize_samples"
# rule_name = "initialize_run"


rule_name = "initialize_run"
rule species_checker:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        component,
        init_complete = component + "/initialize_components_complete",
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
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/initialize_run.py")


# rule get_git_hash_of_serumqc:
#     input:
#         run_info_yaml_path = "run.yaml"
#     output:
#         git_hash = "serumqc/git_hash.txt"
#     run:
#         run_info = datahandling.load_run(input.run_info_yaml_path)
#         shell("git --git-dir {workflow.basedir}/.git rev-parse snakemake 1> {output}")
#         with open(output.git_hash, "r") as git_info:
#             git_hash = git_info.readlines()[0].strip()
#         run_info["run"]["git_hash"] = git_hash
#         datahandling.save_run(run_info, input.run_info_yaml_path)


# rule get_conda_env:
#     input:
#         git_hash = "serumqc/git_hash.txt",
#         run_info_yaml_path = "run.yaml"
#     output:
#         conda_yaml = "serumqc/conda.yaml"
#     run:
#         run_info = datahandling.load_run(input.run_info_yaml_path)
#         shell("conda env export 1> {output}")
#         run_info["run"]["conda_env"] = datahandling.load_yaml(output.conda_yaml)
#         datahandling.save_run(run_info, input.run_info_yaml_path)


# rule add_components_data_entry:
#     input:
#         git_hash = "serumqc/git_hash.txt",
#         conda_yaml = "serumqc/conda.yaml",
#     output:
#         components_db = ""

rule create_end_file:
    input:
        "run.yaml"
    output:
        rules.all.input
    shell:
        """
        bash run_cmd_serumqc.sh
        touch {output}
        """

# can break this down to 2 parts where you create the sample_sheet in one and then prep for run with the other
# rule start_run:
#     input:
#         "run.yaml"
#     output:
#         touch("run_started")
#     shell:
#         "bash run_cmd_serumqc.sh"