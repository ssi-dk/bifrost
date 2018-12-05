import os
import sys
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "../scripts"))
import datahandling


configfile: "../run_config.yaml"
# requires --config R1_reads={read_location},R2_reads={read_location}
sample = config["Sample"]
global_threads = config["threads"]
global_memory_in_GB = config["memory"]

db_sample = datahandling.load_sample(sample)

R1 = db_sample["reads"]["R1"]
R2 = db_sample["reads"]["R2"]

component = "resfinder_db"
version = "1.0"

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


rule_name = "ariba_resfinder"
rule ariba_resfinder:
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
        reads = (R1, R2)
    output:
        folder = directory(rules.setup.params.folder + "/ariba_resfinder")
    params:
        database = os.path.join(os.path.dirname(workflow.snakefile), config["ariba_resfinder_database"])
    conda:
        "../envs/ariba.yaml"
    shell:
        "ariba run {params.database} {input.reads[0]} {input.reads[1]} {output.folder} --tmp_dir /scratch > {log.out_file} 2> {log.err_file}"


rule_name = "datadump_resfinder_db"
rule datadump_analysis:lo
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
        rules.ariba_resfinder.output.folder,
    output:
        summary = touch(rules.all.input)
    params:
        folder = rules.setup.params.folder,
        sample = db_sample.get("name", "ERROR") + "__" + component + ".yaml",
        version = version
    conda:
        "../envs/python_packages.yaml"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_resfinder_db.py")
