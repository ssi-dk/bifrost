import os
import sys
import traceback
import shutil
from bifrostlib import datahandling
from bifrostlib import check_requirements

component = "ariba_resfinder"

configfile: "../config.yaml"  # Relative to run directory
num_of_threads, memory_in_GB = config["threads"], config["memory"]
sample = config["Sample"]

sample_file = sample
db_sample = datahandling.load_sample(sample_file)

component_file = "../components/" + component + ".yaml"
if not os.path.isfile(component_file):
    shutil.copyfile(os.path.join(os.path.dirname(workflow.snakefile), "config.yaml"), component_file)
db_component = datahandling.load_component(component_file)
singularity: db_component["dockerfile"]

sample_component_file = db_sample["name"] + "__" + component + ".yaml"
db_sample_component = datahandling.load_sample_component(sample_component_file)

if "reads" in db_sample:
    reads = R1, R2 = db_sample["reads"]["R1"], db_sample["reads"]["R2"]
else:
    reads = R1, R2 = ("/dev/null", "/dev/null")

onsuccess:
    print("Workflow complete")
    datahandling.update_sample_component_success(sample_component_file, component)


onerror:
    print("Workflow error")
    datahandling.update_sample_component_failure(sample_component_file, component)


rule all:
    input:
        component + "/" + component + "_complete"


rule setup:
    output:
        init_file = touch(temp(component + "/" + component + "_initialized")),
    params:
        folder = component


rule_name = "check_requirements"
rule check_requirements:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = rules.setup.output.init_file,
    output:
        check_file = rules.setup.params.folder + "/requirements_met",
    params:
        sample_file = sample_file,
        component_file = component_file,
        sample_component_file = sample_component_file
    run:
        check_requirements.script__initialization(params.sample_file, params.component_file, params.sample_component_file, output.check_file, log.out_file, log.err_file)


rule_name = "ariba_resfinder"
rule ariba_resfinder:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        check_file = rules.check_requirements.output.check_file,
        folder = rules.setup.output.init_file,
        reads = (R1, R2)
    output:
        complete = rules.setup.params.folder + "/resistance/report.tsv"
    params:
        folder = component,
        sample_file = sample_file,
        component_file = component_file
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/run_ariba_resfinder.py")


rule_name = "datadump"
rule datadump:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        rules.ariba_resfinder.output.complete,
    output:
        complete = rules.all.input
    params:
        folder = rules.setup.params.folder,
        sample_file = sample_file,
        component_file = component_file,
        sample_component_file = sample_component_file
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
