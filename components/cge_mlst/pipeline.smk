import os
import sys
import traceback
import shutil
from bifrostlib import datahandling
from bifrostlib import check_requirements

component = "cge_mlst"

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
        requirements_file = component_file_name
    output:
        check_file = rules.setup.params.folder + "/requirements_met",
    params:
        component = component_file_name,
        sample = sample,
        sample_component = sample_component_file_name
    run:
        check_requirements.script__initialization(input.requirements_file, params.component, params.sample, params.sample_component, output, log.out_file, log.err_file)


rule_name = "cge_mlst"
rule cge_mlst:
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
        folder = rules.setup.output.init_file,
        reads = (R1, R2)
    output:
        complete = rules.setup.params.folder + "/mlst_complete"
    params:
        sample = sample
    run:
        try:
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            species = db_sample["properties"]["species"]

            mlst_species = []
            if species not in db_component["mlst_species_mapping"]:
                datahandling.log(log_out, "cge mlst species: {}\n".format(mlst_species_DB_name))
                shell("touch no_mlst_species_DB")
            else:
                mlst_database_path = os.path.join(os.path.dirname(workflow.snakefile), db_component["mlst_database_path"])
                mlst_species = db_component["mlst_species_mapping"][species]
                for mlst_entry in mlst_species:
                    datahandling.log(log_out, "mlst {} on species: {}\n".format(mlst_entry, species))
                    shell("mkdir {}; mlst.py -x -matrix -s {} -p {} -mp kma -i {} {} -o {} 1> {} 2> {}".format(mlst_entry, mlst_entry, mlst_database_path, input.reads[0], input.reads[1], mlst_entry, log.out_file, log.err_file))
            shell("touch {}".format(output.complete))
            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(traceback.format_exc()))

rule_name = "datadump"
rule datadumpt:
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
        rules.cge_mlst.output.complete,
    output:
        summary = touch(rules.all.input)
    params:
        folder = rules.setup.params.folder,
        sample = db_sample.get("name", "ERROR") + "__" + component + ".yaml",
        sample_yaml = sample
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
