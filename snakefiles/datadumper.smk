import os
import sys
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "../scripts"))
import datahandling

configfile: os.path.join(os.path.dirname(workflow.snakefile), "../config.yaml")
# requires --config R1_reads={read_location},R2_reads={read_location}
sample = config["Sample"]

global_threads = config["threads"]
global_memory_in_GB = config["memory"]

config_sample = datahandling.load_sample(sample)

R1 = config_sample["sample"]["R1"]
R2 = config_sample["sample"]["R2"]

# my understanding is all helps specify final output
component = "datadumper"

onsuccess:
    print("Workflow complete")
    datahandling.update_sample_component_success(sample + "__" + component + ".yaml")

onerror:
    print("Workflow error")
    datahandling.update_sample_component_failure(sample + "__" + component + ".yaml")


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
        directory = directory("datadumper")
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
