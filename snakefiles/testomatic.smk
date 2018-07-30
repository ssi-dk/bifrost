import os
import sys
import datahandling

configfile: "../serumqc_config.yaml"
sample = config["Sample"]
global_threads = config["threads"]
global_memory_in_GB = config["memory"]

sample_doc = datahandling.load_sample(sample)

component = "testomatic"


onsuccess:
    print("Workflow complete")
    datahandling.update_sample_component_success(
        config_sample.get("name", "ERROR") + "__" + component + ".yaml")


onerror:
    print("Workflow error")
    datahandling.update_sample_component_failure(
        config_sample.get("name", "ERROR") + "__" + component + ".yaml")


rule all:
    input:
        component + "/" + component + "_complete"


rule setup:
    output:
        folder = directory(component)
    shell:
        "mkdir {output}"


rule_name = "test_testomatic"
rule test_testomatic:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = rules.setup.output.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.output.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.output.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        qcquickie = "qcquickie.qcquickie_complete",  # Depends on qcquickie
        folder = rules.setup.output,
    output:
        test_results = rules.setup.output.folder + "/test_results.yaml",
    params:
        sample = sample
    script:
        os.path.join(os.path.dirname(workflow.snakefile),
                     "../scripts/testomatic.py")

rule_name = "datadump_testomatic"
rule datadump_testomatic:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = rules.setup.output.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.output.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.output.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = rules.test_testomatic.output,
    output:
        summary = touch(rules.all.input)
    params:
        sample = config_sample.get("name", "ERROR") + \
            "__" + component + ".yaml",
    script:
        os.path.join(os.path.dirname(workflow.snakefile),
                     "../scripts/datadump_testomatic.py")
