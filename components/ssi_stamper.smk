import os
import sys
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "../scripts"))
import datahandling

configfile: "../run_config.yaml"
# requires --config R1_reads={read_location},R2_reads={read_location}

global_threads = config["threads"]
global_memory_in_GB = config["memory"]

sample = config["Sample"]
component = "ssi_stamper"
config_sample = datahandling.load_sample(sample)
sample_component = config_sample["name"] + "__" + component + ".yaml"

R1 = config_sample["reads"]["R1"]
R2 = config_sample["reads"]["R2"]


onsuccess:
    print("Workflow complete")
    datahandling.update_sample_component_success(config_sample.get("name", "ERROR") + "__" + component + ".yaml", component)


onerror:
    print("Workflow error")
    datahandling.update_sample_component_failure(config_sample.get("name", "ERROR") + "__" + component + ".yaml", component)


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
        requirements_file = os.path.join(os.path.dirname(workflow.snakefile), component + ".yaml")
    output:
        check_file = rules.setup.params.folder + "/requirements_met",
    params:
        sample = sample,
        sample_component = sample_component
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/check_requirements.py")


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
        sample = config_sample,
        sample_yaml = sample,
        sample_component = config_sample.get("name", "ERROR") + \
            "__" + component + ".yaml"
    input:
        check_file = rules.check_requirements.output.check_file,
    output:
        complete = touch(rules.all.input)
    script:
        os.path.join(os.path.dirname(workflow.snakefile),
                     "../scripts/ssi_stamper.py")
