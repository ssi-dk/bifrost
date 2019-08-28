#- Templated section: start ------------------------------------------------------------------------
import os
import sys
import traceback
import shutil
from bifrostlib import datahandling
from bifrostlib import check_requirements

configfile: "../config.yaml"  # Relative to run directory
threads, memory_in_GB = config["threads"], config["memory"]
sample = config["Sample"]

sample_file = sample
#db_sample = datahandling.load_sample(sample_file)
db_sample = datahandling.get_sample(sample_id=config["sample_id"])
db_component = datahandling.get_component(component_id=config["component_id"])
db_sample_component = datahandling.get_sample_component(sample_id=config["sample_id"], component_id=config["component_id"])

# TODO: Update db code to be ID based, component has no ID to pass right now and name/version are used instead
# component_file = os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")
# component_config = datahandling.load_yaml(component_file)
# db_component = datahandling.get_components(component_names=[component_config["name"]], component_versions=[component_config["version"]])[0]

singularity: db_component["dockerfile"]

sample_component_file = db_sample["name"] + "__" + db_component["name"] + ".yaml"
db_sample_component = datahandling.load_sample_component(sample_component_file)

if "reads" in db_sample:
    reads = R1, R2 = db_sample["reads"]["R1"], db_sample["reads"]["R2"]
else:
    reads = R1, R2 = ("/dev/null", "/dev/null")

onsuccess:
    print("Workflow complete")
    datahandling.update_sample_component_success(sample_component_file, db_component["name"])


onerror:
    print("Workflow error")
    datahandling.update_sample_component_failure(sample_component_file, db_component["name"])


rule all:
    input:
        db_component["name"] + "/component_complete"


rule setup:
    output:
        init_file = touch(temp(db_component["name"] + "/" + db_component["name"] + "_initialized")),
    params:
        folder = db_component["name"]


rule_name = "check_requirements"
rule check_requirements:
    message:
        "Running step:" + rule_name
    threads:
        threads
    resources:
        memory_in_GB
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
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
#- Templated section: end --------------------------------------------------------------------------

#* Dynamic section: end ****************************************************************************
rule_name = "contaminant_check__classify_reads_kraken_minikraken_db"
rule contaminant_check__classify_reads_kraken_minikraken_db:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        threads
    resources:
        memory_in_GB
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
        threads
    resources:
        memory_in_GB
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
#* Dynamic section: end ****************************************************************************

#- Templated section: start ------------------------------------------------------------------------
rule_name = "datadump"
rule datadump:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        threads
    resources:
        memory_in_GB
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    input:
        #* Dynamic section: start ******************************************************************
        rules.contaminant_check__determine_species_bracken_on_minikraken_results.output.bracken  # Needs to be output of final rule
        #* Dynamic section: end ********************************************************************
    output:
        complete = rules.all.input
    params:
        folder = rules.setup.params.folder,
        sample_file = sample_file,
        component_file = component_file,
        sample_component_file = sample_component_file
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
#- Templated section: end --------------------------------------------------------------------------
