#- Templated section: start ------------------------------------------------------------------------
import os
import sys
import traceback
import shutil
from bifrostlib import datahandling

configfile: "../config.yaml"  # Relative to run directory
num_of_threads, memory_in_GB = config["threads"], config["memory"]

bifrost_sampleComponentObj = datahandling.SampleComponentObj()
(sample_db, component_db) = bifrost_sampleComponentObj.load(config["sample_id"], config["component_id"])

singularity: component_db["dockerfile"]

reads = bifrost_sampleComponentObj.get_reads()

onsuccess:
    bifrost_sampleComponentObj.success()


onerror:
    bifrost_sampleComponentObj.failure()


rule all:
    input:
        component_db["name"] + "/datadump_complete"  # file is defined by datadump function


rule setup:
    output:
        init_file = touch(temp(component_db["name"] + "/" + component_db["name"] + "_initialized")),
    params:
        folder = component_db["name"]


rule_name = "check_requirements"
rule check_requirements:
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = component_db["name"] + "/log/" + rule_name + ".out.log",
        err_file = component_db["name"] + "/log/" + rule_name + ".err.log",
    benchmark:
        component_db["name"] + "/benchmarks/" + rule_name + ".benchmark"
    input:
        folder = rules.setup.output.init_file,
    output:
        check_file = component_db["name"] + "/requirements_met",
    params:
        bifrost_sampleComponentObj
    run:
        bifrost_sampleComponentObj.check_requirements()
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
        out_file = component_db["name"] + "/log/" + rule_name + ".out.log",
        err_file = component_db["name"] + "/log/" + rule_name + ".err.log",
    benchmark:
        component_db["name"] + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        rules.check_requirements.output.check_file,
        reads
    output:
        stats_file = component_db["name"] + "/stats.txt"
    params:
        adapters = component_db["resources"]["adapters_fasta"]  # This is now done to the root of the continuum container
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
        out_file = component_db["name"] + "/log/" + rule_name + ".out.log",
        err_file = component_db["name"] + "/log/" + rule_name + ".err.log",
    benchmark:
        component_db["name"] + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        stats_file = rules.setup__filter_reads_with_bbduk.output.stats_file,
    params:
        sampleComponentObj=bifrost_sampleComponentObj
    output:
        file = component_db["name"] + "/has_min_num_of_reads"
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
        out_file = component_db["name"] + "/log/" + rule_name + ".out.log",
        err_file = component_db["name"] + "/log/" + rule_name + ".err.log",
    benchmark:
        component_db["name"] + "/benchmarks/" + rule_name + ".benchmark"
    input:
        #* Dynamic section: start ******************************************************************
        rules.greater_than_min_reads_check.output.file  # Needs to be output of final rule
        #* Dynamic section: end ********************************************************************
    output:
        complete = rules.all.input
    params:
        sampleComponentObj=bifrost_sampleComponentObj
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
#- Templated section: end --------------------------------------------------------------------------
