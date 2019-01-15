import os
import sys
from bifrostlib import datahandling

component = "analyzer"  # Depends on component name, should be same as folder

configfile: "../run_config.yaml"  # Relative to run directory
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


rule_name = "ariba_resfinder"
rule ariba_resfinder:
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
        folder = directory(rules.setup.params.folder + "/ariba_resfinder")
    params:
        database = os.path.join(os.path.dirname(workflow.snakefile), db_component["ariba_resfinder_database"])
    shell:
        "ariba run {params.database} {input.reads[0]} {input.reads[1]} {output.folder} --tmp_dir /scratch > {log.out_file} 2> {log.err_file}"


rule_name = "abricate_on_ariba_resfinder"
rule abricate_on_ariba_resfinder:
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
        contigs = rules.ariba_resfinder.output.folder
    output:
        report = rules.setup.params.folder + "/abricate_on_resfinder_from_ariba.tsv",
    params:
        database = os.path.join(os.path.dirname(workflow.snakefile), db_component["abricate_resfinder_database"])
    shell:
        """
        if [[ -e {input.contigs}/assemblies.fa.gz ]] && [[ -n $(gzip -cd {input.contigs}/assemblies.fa.gz | head -c1) ]];
        then abricate --datadir {params.database} --db . {input.contigs}/assemblies.fa.gz > {output.report} 2> {log.err_file};
        else touch {output.report};
        fi;
        """


rule_name = "ariba_plasmidfinder"
rule ariba_plasmidfinder:
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
        folder = directory(rules.setup.params.folder + "/ariba_plasmidfinder")
    params:
        database = os.path.join(os.path.dirname(workflow.snakefile), db_component["ariba_plasmidfinder_database"])
    shell:
        "ariba run {params.database} {input.reads[0]} {input.reads[1]} {output.folder} --tmp_dir /scratch > {log.out_file} 2> {log.err_file}"


rule_name = "abricate_on_ariba_plasmidfinder"
rule abricate_on_ariba_plasmidfinder:
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
        folder = rules.ariba_plasmidfinder.output
    output:
        report = rules.setup.params.folder + "/abricate_on_plasmidfinder_from_ariba.tsv",
    params:
        database = os.path.join(os.path.dirname(workflow.snakefile), db_component["abricate_plasmidfinder_database"])
    shell:
        """
        if [[ -e {input.folder}/assemblies.fa.gz ]] && [[ -n $(gzip -cd {input.folder}/assemblies.fa.gz | head -c1) ]];
        then abricate --datadir {params.database} --db . {input.folder}/assemblies.fa.gz > {output.report} 2> {log.err_file};
        else touch {output.report};
        fi;
        """


rule_name = "ariba_mlst"
rule ariba_mlst:
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
        folder = directory(rules.setup.params.folder + "/ariba_mlst")
    params:
        sample = sample
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
            datahandling.log(log_err, str(traceback.format_exc()))

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
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
