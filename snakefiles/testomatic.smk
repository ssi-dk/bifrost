import os
import sys
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "../scripts"))
import datahandling

configfile: "../serumqc_config.yaml"
sample = config["Sample"]
global_threads = config["threads"]
global_memory_in_GB = config["memory"]

config_sample = datahandling.load_sample(sample)

component = "testomatic"


onsuccess:
    print("Workflow complete")
    datahandling.update_sample_component_success(
        config_sample.get("name", "ERROR") + "__" + component + ".yaml")


onerror:
    print("Workflow error")
    datahandling.update_sample_component_failure(
        config_sample.get("name", "ERROR") + "__" + component + ".yaml")

# ruleorder: setup > test_testomatic > datadump_testomatic > all

rule all:
    input:
        component + "/" + component + "_complete"


rule setup:
    output:
        init_file = touch(temp(component + "/" + component + "_initialized")),
    params:
        folder = component,
    # shell:
    #     """
    #     mkdir {output}
    #     """


rule_name = "run_testomatic"
rule run_testomatic:
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
        sample_component = config_sample.get("name", "ERROR") + \
            "__" + component + ".yaml",
        assemblatron = config_sample["name"] + "__assemblatron.yaml", # optional
    input:
        qcquickie = "qcquickie/qcquickie_complete",  # Depends on qcquickie
        qcquickie_yaml = config_sample["name"] + "__qcquickie.yaml",
        folder = rules.setup.output,
    output:
        test_results = rules.setup.params.folder + "/test_results.yaml",
        complete = rules.all.input
    script:
        os.path.join(os.path.dirname(workflow.snakefile),
                     "../scripts/testomatic.py")
