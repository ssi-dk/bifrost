#- Templated section: start ------------------------------------------------------------------------
import os
import sys
import traceback
import shutil
from bifrostlib import datahandling
from bifrostlib import check_requirements

configfile: "../config.yaml"  # Relative to run directory
num_of_threads, memory_in_GB = config["threads"], config["memory"]
sample = config["Sample"]

sample_file = sample
sample_db = datahandling.load_sample(sample_file)

component_file = os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")
component_db = datahandling.load_component(component_file)

singularity: component_db["dockerfile"]

sample_component_file = sample_db["name"] + "__" + component_db["name"] + ".yaml"
db_sample_component = datahandling.load_sample_component(sample_component_file)

if "reads" in sample_db:
    reads = R1, R2 = sample_db["reads"]["R1"], sample_db["reads"]["R2"]
else:
    reads = R1, R2 = ("/dev/null", "/dev/null")

onsuccess:
    print("Workflow complete")
    datahandling.update_sample_component_success(sample_component_file, component_db["name"])


onerror:
    print("Workflow error")
    datahandling.update_sample_component_failure(sample_component_file, component_db["name"])


rule all:
    input:
        component_db["name"] + "/" + component_db["name"] + "_complete"


rule setup:
    output:
        init_file = touch(temp(component_db["name"] + "/" + component_db["name"] + "_initialized")),
    params:
        folder = component_db["name"]


rule_name = "check_requirements"
rule check_requirements:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
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
#- Templated section: end --------------------------------------------------------------------------

#* Dynamic section: start **************************************************************************
rule_name = "rule_1"
rule rule_1:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
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
        file = rules.setup.params.folder + "/data.json"
    params:
        folder = component_db["name"],
        sample_file = sample_file,
        component_file = component_file
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/run_rule_name.py")
#* Dynamic section: end ****************************************************************************

#- Templated section: start ------------------------------------------------------------------------
rule_name = "datadump"
rule datadump:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        #* Dynamic section: start ******************************************************************
        rules.rule_1.output.complete  # Needs to be output of final rule
        #* Dynamic section: end ********************************************************************
    output:
        complete = rules.all.input
    params:
        folder = rules.setup.params.folder,
        sample_file = sample_file,
        component_file = component_file,
        sample_component_file = sample_component_file
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
#- Templated section: end --------------------------------------------------------------------------
