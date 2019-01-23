import os
import sys
import Bio.SeqIO
import traceback
import shutil
import pandas
from bifrostlib import datahandling

component = "whats_my_species"  # Depends on component name, should be same as folder

configfile: "../config.yaml"  # Relative to run directory
global_threads = config["threads"]
global_memory_in_GB = config["memory"]
sample = config["Sample"]

sample_file_name = sample
db_sample = datahandling.load_sample(sample_file_name)

component_file_name = "../components/" + component + ".yaml"
if not os.path.isfile(component_file_name):
    shutil.copyfile(os.path.join(os.path.dirname(workflow.snakefile), "config.yaml"), component_file_name)
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
        requirements_file = component_file_name
    output:
        check_file = rules.setup.params.folder + "/requirements_met",
    params:
        component = component_file_name,
        sample = sample,
        sample_component = sample_component_file_name
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../common/check_requirements.py")


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
        db = os.path.join(os.path.dirname(workflow.snakefile), db_component["kraken_database"])
    shell:
        "kraken --threads {threads} -db {params.db} --fastq-input {input.reads} 2> {log.err_file} | kraken-report -db {params.db} 1> {output.kraken_report}"


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
        kmer_dist = os.path.join(os.path.dirname(workflow.snakefile), db_component["kraken_kmer_dist"])
    shell:
        """
        est_abundance.py -i {input.kraken_report} -k {params.kmer_dist} -o {output.bracken} 1> {log.out_file} 2> {log.err_file}
        sort -r -t$'\t' -k7 {output.bracken} -o {output.bracken}
        """

rule_name = "species_check__set_species"
rule species_check__set_species:
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
        bracken = rules.contaminant_check__determine_species_bracken_on_minikraken_results.output.bracken,
    output:
        species = rules.setup.params.folder + "/species.txt",
    params:
        sample = sample,
    run:
        try:
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            sample_db = datahandling.load_sample(sample)
            with open(output.species, "w") as species_file:
                df = pandas.read_table(input.bracken)
                # This try-except will avoid a crash when no species is found by bracken.
                try:
                    sample_db["properties"]["detected_species"] = df["name"].iloc[0]
                except IndexError:
                    sample_db["properties"]["detected_species"] = None
                sample_db["properties"]["provided_species"] = sample_db["properties"].get("provided_species",)
                if sample_db["properties"]["provided_species"] is not None:
                    sample_db["properties"]["species"] = sample_db["properties"]["provided_species"]
                else:
                    sample_db["properties"]["species"] = sample_db["properties"]["detected_species"]
            datahandling.save_sample(sample_db, sample)
        except Exception as e:
            datahandling.log(log_err, str(traceback.format_exc()))

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
        species = rules.species_check__set_species.output.species,
    output:
        summary = touch(rules.all.input)
    params:
        sample = db_sample.get("name", "ERROR") + "__" + component + ".yaml",
        folder = rules.setup.params.folder
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
