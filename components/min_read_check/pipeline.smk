#- Templated section: start ------------------------------------------------------------------------
import os
import sys
import traceback
import shutil
from bifrostlib import datahandling
from bifrostlib import check_requirements

configfile: "../config.yaml"  # Relative to run directory
num_of_threads, memory_in_GB = config["threads"], config["memory"]
memory_in_GB = config["memory"]
sample = config["Sample"]

sample_file = sample

# TODO: Update db code to be ID based, component has no ID to pass right now and name/version are used instead. Probably a helper function to load the 3
db_sample = datahandling.load_sample(sample_file)
component_file = os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")
component_config = datahandling.load_yaml(component_file)
db_component = datahandling.get_components(component_names=[component_config["name"]], component_versions=[component_config["version"]])[0]
sample_component_file = db_sample["name"] + "__" + db_component["name"] + ".yaml"

singularity: db_component["dockerfile"]

if "reads" in db_sample:
    reads = R1, R2 = db_sample["reads"]["R1"], db_sample["reads"]["R2"]
else:
    reads = R1, R2 = ("/dev/null", "/dev/null")


onsuccess:
    print("Workflow complete")
    datahandling.update_sample_component_success(sample_component_file)


onerror:
    print("Workflow error")
    datahandling.update_sample_component_failure(sample_component_file)


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
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = db_component["name"] + "/log/" + rule_name + ".out.log",
        err_file = db_component["name"] + "/log/" + rule_name + ".err.log",
    benchmark:
        db_component["name"] + "/benchmarks/" + rule_name + ".benchmark"
    input:
        folder = rules.setup.output.init_file,
    output:
        check_file = db_component["name"] + "/requirements_met",
    params:
        sample_file = sample_file,
        component_file = component_file,
        sample_component_file = sample_component_file
    run:
        check_requirements.script__initialization(params.sample_file, params.component_file, params.sample_component_file, output.check_file, log.out_file, log.err_file)
#- Templated section: end --------------------------------------------------------------------------

#* Dynamic section: end ****************************************************************************
rule_name = "setup__filter_reads_with_bbduk"
rule setup__filter_reads_with_bbduk:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = db_component["name"] + "/log/" + rule_name + ".out.log",
        err_file = db_component["name"] + "/log/" + rule_name + ".err.log",
    benchmark:
        db_component["name"] + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        rules.check_requirements.output.check_file,
        reads = (R1, R2)
    output:
        stats_file = db_component["name"] + "/stats.txt"
    params:
        adapters = db_component["resources"]["adapters_fasta"]  # This is now done to the root of the continuum container
    shell:
        "bbduk.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.reads[0]} in2={input.reads[1]} ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tbo qtrim=r minlength=30 1> {log.out_file} 2> {output.stats_file}"


rule_name = "greater_than_min_reads_check"
rule greater_than_min_reads_check:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = db_component["name"] + "/log/" + rule_name + ".out.log",
        err_file = db_component["name"] + "/log/" + rule_name + ".err.log",
    benchmark:
        db_component["name"] + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        stats_file = rules.setup__filter_reads_with_bbduk.output.stats_file,
    params:
        sample_file = sample_file,
        component_file = component_file
    output:
        file = db_component["name"] + "/has_min_num_of_reads"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/rule__greater_than_min_reads_check.py")
#* Dynamic section: end ****************************************************************************

#- Templated section: start ------------------------------------------------------------------------
rule_name = "datadump"
rule datadump:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = db_component["name"] + "/log/" + rule_name + ".out.log",
        err_file = db_component["name"] + "/log/" + rule_name + ".err.log",
    benchmark:
        db_component["name"] + "/benchmarks/" + rule_name + ".benchmark"
    input:
        #* Dynamic section: start ******************************************************************
        rules.greater_than_min_reads_check.output.file  # Needs to be output of final rule
        #* Dynamic section: end ********************************************************************
    output:
        complete = rules.all.input
    params:
        folder = db_component["name"],
        sample_file = sample_file,
        component_file = component_file,
        sample_component_file = sample_component_file
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
#- Templated section: end --------------------------------------------------------------------------
