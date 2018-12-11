import os
import sys
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "../../scripts"))
import datahandling from bifrostlib


component = "sp_cdiff_fbi"  # Depends on component name, should be same as folder

configfile: "../run_config.yaml"  # Relative to run directory
global_threads = config["threads"]
global_memory_in_GB = config["memory"]
sample = config["Sample"]

sample_file_name = sample
db_sample = datahandling.load_sample(sample_file_name)

component_file_name = os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")
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
        requirements_file = os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")
    output:
        check_file = rules.setup.params.folder + "/requirements_met",
    params:
        sample = sample,
        sample_component = sample_component_file_name
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../../scripts/check_requirements.py")



rule_name = "cdiff_analysis"
rule cdiff_analysis:
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
        rules.check_requirements.output.check_file,
    output:
        folder = directory(rules.setup.params.folder + "/cdiff_analysis")
    params:
        sample = sample
    shell:
        # This and the name are the only things that need to be changed, please ask for assistance
        "/srv/data/tools/git.repositories/SSI-scripts/qcscripts/qcfindgene.sh"


"""
This should be done but is being hacked around for now
"""
# rule_name = "datadump_analyzer"
# rule datadump_analysis:
#     # Static
#     message:
#         "Running step:" + rule_name
#     threads:
#         global_threads
#     resources:
#         memory_in_GB = global_memory_in_GB
#     log:
#         out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
#         err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
#     benchmark:
#         rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
#     # Dynamic
#     input:
#         rules.cdiff_analysis.output.folder,
#     output:
#         summary = touch(rules.all.input)
#     params:
#         folder = rules.setup.params.folder,
#         sample = db_sample.get("name", "ERROR") + "__" + component + ".yaml",
#     conda:
#         "../envs/python_packages.yaml"
#     shell:
#         os.path.join(os.path.dirname(workflow.snakefile), "../../scripts/datadump_analyzer.py")

rule temp_finish:
    input:
        rules.cdiff_analysis.output.folder
    output:
        summary = touch(rules.all.input)
    shell:
        "touch {output.summary}"