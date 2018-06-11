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


onsuccess:
    print("Workflow complete")
    shell("mv datadumper/summary.yaml sample.yaml")
    with open(sample, "r") as sample_yaml:
        config_sample = yaml.load(sample_yaml)
    if "datadumper" not in config_sample["sample"]["components"]["success"]:
        config_sample["sample"]["components"]["success"].append("datadumper")
    if "datadumper" in config_sample["sample"]["components"]["failure"]:
        config_sample["sample"]["components"]["failure"].remove("datadumper")
    with open(sample, "w") as output_file:
        yaml.dump(config_sample, output_file)

onerror:
    print("Workflow error")
    if "datadumper" in config_sample["sample"]["components"]["success"]:
        config_sample["sample"]["components"]["success"].remove("datadumper")
    if "datadumper" not in config_sample["sample"]["components"]["failure"]:
        config_sample["sample"]["components"]["failure"].append("datadumper")
    with open(sample, "w") as output_file:
        yaml.dump(config_sample, output_file)


rule all:
    input:
        "datadumper",
        "datadumper/summary.yaml"


rule setup:
    message:
        "Running step: {rule}"
    output:
        dir = "datadumper"
    shell:
        "mkdir {output}"


rule datadump_qcquickie:
    message:
        "Running step: {rule}"
    input:
        datadumper = "datadumper",
        folder = "qcquickie",
    output:
        summary = "datadumper/qcquickie.yaml"
    params:
        sample = config_sample
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    # conda:
    #     "../envs/fastqc.yaml"
    group:
        "datadumper"
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
        folder = "assembly",
    output:
        summary = "datadumper/assembly.yaml"
    params:
        sample = config_sample
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    # conda:
    #     "../envs/fastqc.yaml"
    group:
        "datadumper"
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
        folder = "analysis",
    output:
        summary = "datadumper/analysis.yaml"
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    # conda:
    #     "../envs/fastqc.yaml"
    group:
        "datadumper"
    log:
        "datadumper/log/datadump_analysis.log"
    benchmark:
        "datadumper/benchmarks/datadump_analysis.benchmark"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_analysis.py")


rule combine_datadumps:
    message:
        "Running step: {rule}"
    input:
        qcquickie_summary = "datadumper/qcquickie.yaml",
        assembly_summary = "datadumper/assembly.yaml",
    output:
        summary = "datadumper/summary.yaml",
    params:
        sample_yaml = "sample.yaml",
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    # conda:
    #     "../envs/fastqc.yaml"
    group:
        "datadumper"
    log:
        "datadumper/log/combine_datadumps.log"
    benchmark:
        "datadumper/benchmarks/combine_datadumps.benchmark"
    # shell:
    #     "cat {input.qcquickie_summary} {input.assembly_summary} {params.sample_yaml} > {output.summary}"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_combine.py")
