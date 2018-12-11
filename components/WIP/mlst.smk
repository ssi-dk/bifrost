import os
import sys
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "../scripts"))
import datahandling from bifrostlib


configfile: "../run_config.yaml"
# requires --config R1_reads={read_location},R2_reads={read_location}
sample = config["Sample"]
global_threads = config["threads"]
global_memory_in_GB = config["memory"]

db_sample = datahandling.load_sample(sample)

R1 = db_sample["reads"]["R1"]
R2 = db_sample["reads"]["R2"]

component = ""


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


rule_name = "check_required_components"
rule check_required_components:
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
        check_file = rules.setup.params.folder + "/required_components_present",
    run:
        try:
            check_file = str(output.check_file)
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            required_components = ["whats_my_species"]
            required_components_success = True
            for component in required_components:
                if not datahandling.sample_component_success(db_sample.get("name", "ERROR") + "__" + component + ".yaml", component):
                    required_components_successful = False
                    datahandling.log(log_out, "Missing component: {}\n".format(component))
            if required_components_success:
                with open(check_file, "w") as out_file:
                    datahandling.log(log_out, "Required components found: {}\n".format(",".join(required_components)))

            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(e))


rule_name = "ariba_mlst"
rule ariba_mlst:
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
        rules.check_required_components.output.check_file,
        folder = rules.setup.output.init_file,
        reads = (R1, R2)
    output:
        folder = directory(rules.setup.params.folder + "/ariba_mlst")
    params:
        sample = sample
    conda:
        "../envs/ariba.yaml"
    run:
        try:
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            mlst_species_DB = datahandling.get_mlst_species_DB(sample)
            datahandling.log(log_out, "mlst species: {}\n".format(mlst_species_DB))
            if mlst_species_DB is None:
                shell("mkdir {}".format(output.folder))
                shell("touch {}/no_mlst_species_DB".format(output.folder))
            else:
                shell("ariba run {} {} {} {} 1> {} 2> {}".format(mlst_species_DB, input.reads[0], input.reads[1], output.folder, log.out_file, log.err_file))
            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(e))

rule_name = "datadump_analyzer"
rule datadump_analysis:
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
        rules.abricate_on_ariba_resfinder.output.report,
        rules.abricate_on_ariba_plasmidfinder.output.report,
        folder = rules.ariba_mlst.output.folder,
    output:
        summary = touch(rules.all.input)
    params:
        folder = rules.setup.params.folder,
        sample = db_sample.get("name", "ERROR") + "__" + component + ".yaml",
    conda:
        "../envs/python_packages.yaml"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_analyzer.py")
