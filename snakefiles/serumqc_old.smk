import re
import pandas
from ruamel.yaml import YAML
import sys


configfile: os.path.join(os.path.dirname(workflow.snakefile), "../config.yaml")
# requires --config R1_reads={read_location},R2_reads={read_location}
sample = config["Sample"]

global_threads = config["threads"]
global_memory_in_GB = config["memory"]

yaml = YAML(typ='safe')
yaml.default_flow_style = False

with open(sample, "r") as yaml_stream:
    config_sample = yaml.load(yaml_stream)

R1 = config_sample["sample"]["R1"]
R2 = config_sample["sample"]["R2"]

# my understanding is all helps specify final output
component = "serumqc_old"

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
        component + "/assembly_complete"


rule setup:
    message:
        "Running step: {rule}"
    output:
        directory = component,
    shell:
        "mkdir {output}"

rule run_serumqc_old:
    message:
        "Running step: {rule}"
    input:
        directory = component,
        reads = (R1, R2)
    output:
        complete = component + "/assembly_complete",
        read_folder = temp(component + "/serumqc_read_folder"),
        folder = component + "/serumqc_old"
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        component + "/log/run_serumqc_old.log"
    benchmark:
        component + "/benchmarks/run_serumqc_old.benchmark"
    shell:
        """
        mkdir {output.read_folder}
        ln {input.reads[0]} {output.read_folder}
        ln {input.reads[1]} {output.read_folder}
        source activate env_serumqc
        serumqc.py -i {output.read_folder} -o {output.folder} -run serumqc_wrapper
        touch {output.complete}
        """

