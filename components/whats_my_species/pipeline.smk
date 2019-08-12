import os
import sys
import traceback
import shutil
import pandas
from bifrostlib import datahandling
from bifrostlib import check_requirements

component = "whats_my_species"  # Depends on component name, should be same as folder

configfile: "../config.yaml"  # Relative to run directory
global_threads = config["threads"]
global_memory_in_GB = config["memory"]
sample = config["Sample"]

sample_file = sample
db_sample = datahandling.load_sample(sample_file)

component_file = "../components/" + component + ".yaml"
if not os.path.isfile(component_file):
    shutil.copyfile(os.path.join(os.path.dirname(workflow.snakefile), "config.yaml"), component_file)
db_component = datahandling.load_component(component_file)
singularity: db_component["dockerfile"]

sample_component_file = db_sample["name"] + "__" + component + ".yaml"
db_sample_component = datahandling.load_sample_component(sample_component_file)

if "reads" in db_sample:
    reads = R1, R2 = db_sample["reads"]["R1"], db_sample["reads"]["R2"]
else:
    reads = R1, R2 = ("/dev/null", "/dev/null")

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
    output:
        check_file = rules.setup.params.folder + "/requirements_met",
    params:
        sample_file = sample_file,
        component_file = component_file,
        sample_component_file = sample_component_file
    run:
        check_requirements.script__initialization(params.sample_file, params.component_file, params.sample_component_file, output.check_file, log.out_file, log.err_file)


rule_name = "contaminant_check__classify_reads_kraken_minikraken_db"
rule contaminant_check__classify_reads_kraken_minikraken_db:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        rules.check_requirements.output.check_file,
        reads = (R1, R2)
    output:
        kraken_report = rules.setup.params.folder + "/kraken_report.txt"
    params:
        db = db_component["kraken_database"]
    shell:
        "kraken --threads {threads} -db {params.db} {input.reads} 2> {log.err_file} | kraken-report -db {params.db} 1> {output.kraken_report}"


rule_name = "contaminant_check__determine_species_bracken_on_minikraken_results"
rule contaminant_check__determine_species_bracken_on_minikraken_results:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        kraken_report = rules.contaminant_check__classify_reads_kraken_minikraken_db.output.kraken_report,
    output:
        bracken = rules.setup.params.folder + "/bracken.txt",
        kraken_report_bracken = rules.setup.params.folder + "/kraken_report_bracken.txt"
    params:
        kmer_dist = db_component["kraken_kmer_dist"]
    shell:
        """
        est_abundance.py -i {input.kraken_report} -k {params.kmer_dist} -o {output.bracken} 1> {log.out_file} 2> {log.err_file}
        sort -r -t$'\t' -k7 {output.bracken} -o {output.bracken}
        """


rule_name = "datadump_whats_my_species"
rule datadump_whats_my_species:
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
        species = rules.contaminant_check__determine_species_bracken_on_minikraken_results.output.kraken_report_bracken,
    output:
        summary = touch(rules.all.input)
    params:
        folder = rules.setup.params.folder,
        sample_file = sample_file,
        component_file = component_file,
        sample_component_file = sample_component_file
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
