import os
import sys
from bifrostlib import datahandling

component = "ssi_stamper"  # Depends on component name, should be same as folder

configfile: "../run_config.yaml"  # Relative to run directory
global_threads = config["threads"]
global_memory_in_GB = config["memory"]
sample = config["Sample"]

sample_file_name = sample
db_sample = datahandling.load_sample(sample_file_name)

component_file_name = "../components/" + component + ".yaml"
if not os.path.isfile(component_file_name):
    shutil.copyfile(os.path.join(os.path.dirname(workflow.snakefile), "config.yaml"), component_file_name)
db_component = datahandling.load_component(component_file_name)

sample_component_file_name = db_sample["name"] + "__" + component + ".yaml"
db_sample_component = datahandling.load_sample_component(sample_component_file_name)

reads = R1, R2 = db_sample["reads"]["R1"], db_sample["reads"]["R2"]

onsuccess:
    print("Workflow complete")
    datahandling.update_sample_component_success(db_sample.get("name", "ERROR") + "__" + component + ".yaml", component)


onerror:
    print("Workflow error")
    datahandling.update_sample_component_failure(db_sample.get("name", "ERROR") + "__" + component + ".yaml", component)


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
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = rules.setup.output.init_file,
        requirements_file = component_file_name
    output:
        check_file = rules.setup.params.folder + "/requirements_met",
    params:
        component = component_file_name,
        sample = sample,
        sample_component = sample_component_file_name
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../common/check_requirements.py")


rule_name = "run_ssi_stamper"
rule run_ssi_stamper:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    params:
        sample = db_sample,
        sample_yaml = sample,
        sample_component = db_sample.get("name", "ERROR") + \
            "__" + component + ".yaml"
    input:
        check_file = rules.check_requirements.output.check_file,
    output:
        complete = touch(rules.all.input)
    script:
        # Should be refactored to a datadump
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/ssi_stamper.py")
