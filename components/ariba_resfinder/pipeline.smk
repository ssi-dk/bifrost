#- Templated section: start ------------------------------------------------------------------------
import os
from bifrostlib import datahandling

configfile: "../config.yaml"  # Relative to run directory
num_of_threads, memory_in_GB = config["threads"], config["memory"]

bifrost_sampleComponentObj = datahandling.SampleComponentObj()
(sample_db, component_db) = bifrost_sampleComponentObj.load(
    config["sample_id"], config["component_id"])

singularity: component_db["dockerfile"]

reads = bifrost_sampleComponentObj.get_reads()

onsuccess:
    bifrost_sampleComponentObj.success()


onerror:
    bifrost_sampleComponentObj.failure()


rule all:
    input:
        # file is defined by datadump function
        component_db["name"] + "/datadump_complete"


rule setup:
    output:
        init_file = touch(
            temp(component_db["name"] + "/" + component_db["name"] + "_initialized")),
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

#* Dynamic section: start **************************************************************************
rule_name = "ariba_resfinder"
rule ariba_resfinder:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        check_file = rules.check_requirements.output.check_file,
        folder = rules.setup.output.init_file,
        reads = reads
    output:
        complete = rules.setup.params.folder + "/resistance/report.tsv"
    params:
        sampleComponentObj = bifrost_sampleComponentObj
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/rule__run_ariba_resfinder.py")


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
        rules.ariba_resfinder.output.file  # Needs to be output of final rule
        #* Dynamic section: end ********************************************************************
    output:
        complete = rules.all.input
    params:
        sampleComponentObj = bifrost_sampleComponentObj
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
#- Templated section: end --------------------------------------------------------------------------
