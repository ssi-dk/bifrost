import re
import pandas
from ruamel.yaml import YAML
import sys


configfile: os.path.join(os.path.dirname(workflow.snakefile), "../config.yaml")

sample = config["Sample"]  # expected input

global_threads = config["threads"]
global_memory_in_GB = config["memory"]

yaml = YAML(typ='safe')
yaml.default_flow_style = False

with open(sample, "r") as sample_yaml:
    config_sample = yaml.load(sample_yaml)


component = "datadumper"

onsuccess:
    print("Workflow complete")
    with open(sample, "r") as sample_yaml:
        config_sample = yaml.load(sample_yaml)
    while component in config_sample["sample"]["components"]["failure"]:
        config_sample["sample"]["components"]["failure"].remove(component)
    if component not in config_sample["sample"]["components"]["success"]:
        config_sample["sample"]["components"]["success"].append(component)
    with open(sample, "w") as output_file:
        yaml.dump(config_sample, output_file)

onerror:
    print("Workflow error")
    with open(sample, "r") as sample_yaml:
        config_sample = yaml.load(sample_yaml)
    while component in config_sample["sample"]["components"]["success"]:
        config_sample["sample"]["components"]["failure"].remove(component)
    if component not in config_sample["sample"]["components"]["failure"]:
        config_sample["sample"]["components"]["success"].append(component)
    with open(sample, "w") as output_file:
        yaml.dump(config_sample, output_file)


rule all:
    input:
        "datadumper",
        "datadumper/qcquickie.yaml",
        "datadumper/assembly.yaml",
        "datadumper/analysis.yaml",


rule setup:
    message:
        "Running step: {rule}"
    output:
        directory = "datadumper"
    shell:
        "mkdir {output}"


rule datadump_qcquickie:
    message:
        "Running step: {rule}"
    input:
        datadumper = "datadumper",
    output:
        summary = "datadumper/qcquickie.yaml"
    params:
        sample = sample,
        folder = "qcquickie",
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        "datadumper/log/datadump_qcquickie.log"
    benchmark:
        "datadumper/benchmarks/datadump_qcquickie.benchmark"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_qcquickie.py")


rule datadump_assembly:
    message:
        "Running step: {rule}"
    input:
        datadumper = "datadumper",
    output:
        summary = "datadumper/assembly.yaml"
    params:
        sample = sample,
        folder = "assembly",
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        "datadumper/log/datadump_assembly.log"
    benchmark:
        "datadumper/benchmarks/datadump_assembly.benchmark"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_assembly.py")


rule datadump_analysis:
    message:
        "Running step: {rule}"
    input:
        datadumper = "datadumper",
    output:
        summary = "datadumper/analysis.yaml"
    params:
        sample = sample,
        folder = "analysis",
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        "datadumper/log/datadump_analysis.log"
    benchmark:
        "datadumper/benchmarks/datadump_analysis.benchmark"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_analysis.py")


# rule combine_datadumps:
#     message:
#         "Running step: {rule}"
#     input:
#         datadumper = "datadumper",
#         qcquickie_summary = "datadumper/qcquickie.yaml",
#         assembly_summary = "datadumper/assembly.yaml",
#     output:
#         summary = "datadumper/summary.yaml",
#     params:
#         sample_yaml = "sample.yaml",
#     threads:
#         global_threads
#     resources:
#         memory_in_GB = global_memory_in_GB
#     log:
#         "datadumper/log/combine_datadumps.log"
#     benchmark:
#         "datadumper/benchmarks/combine_datadumps.benchmark"
#     script:
#         os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_combine.py")
