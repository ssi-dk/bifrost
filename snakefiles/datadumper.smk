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
    config_sample["sample"]["components"]["success"].append("datadumper")
    print(config_sample)
    with open(sample, "w") as output_file:
        yaml.dump(config_sample, output_file)

onerror:
    print("Workflow error")
    config_sample["sample"]["components"]["failure"].append("datadumper")
    with open(sample, "w") as output_file:
        yaml.dump(config_sample, output_file)


rule all:
    input:
        "datadumper",
        "datadumper/summary.yaml"


rule setup:
    output:
        dir = "datadumper"
    group:
        "datadumper"
    shell:
        "mkdir {output}"


rule datadump_qcquickie:
    message:
        "Running step: {rule}"
    input:
        folder = "qcquickie",
    output:
        summary = "qcquickie/qcquickie.yaml"
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
        folder = "assembly",
    output:
        summary = "assembly/qcquickie.yaml"
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
        folder = "analysis",
    output:
        summary = "analysis/analysis.yaml"
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
        folder = "qcquickie/qcquickie.yaml",
    output:
        summary = pipe("datadumper/summary.yaml"),
    params:
        sample = "sample.yaml",
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
    shell:
        "cat {input.folder} {params.sample} 1> {summary.yaml} 2> {log}"
