import os
import sys
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "../../scripts"))
import datahandling


configfile: "../run_config.yaml"
# requires --config R1_reads={read_location},R2_reads={read_location}
sample = config["Sample"]
global_threads = config["threads"]
global_memory_in_GB = config["memory"]

config_sample = datahandling.load_sample(sample)

R1 = config_sample["reads"]["R1"]
R2 = config_sample["reads"]["R2"]

component = "sp_cdiff_fbi"


onsuccess:
    print("Workflow complete")
    datahandling.update_sample_component_success(config_sample.get("name", "ERROR") + "__" + component + ".yaml", component)


onerror:
    print("Workflow error")
    datahandling.update_sample_component_failure(config_sample.get("name", "ERROR") + "__" + component + ".yaml", component)


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
    params:
        requirements = [{"test":"value"}]
    output:
        check_file = rules.setup.params.folder + "/requirements_met",
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/check_requirements.py")



rule_name = "cdiff_analysis"
rule cdiff_analysis:
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

    output:
        folder = directory(rules.setup.params.folder + "/ done_qcfindgene")
    params:
        sample = sample
    shell:
        "/srv/data/tools/git.repositories/SSI-scripts/qcscripts/qcfindgene.sh"
        #sample_analyzer_db = datahandling.load_sample_component(sample + "__analyzer.yaml")
        # load sample datahandline.load_sample("../sample.yaml")
    run:
        # try:
        #     config_sample = config_sample
        #     log_out = str(log.out_file)
        #     log_err = str(log.err_file)
            
        #     sample_analyzer_db = datahandling.load_sample_component(sample + "__analyzer.yaml")
        #     mlst_string = sample_analyzer_db["summary"].get("mlst", None)
        #     # "ST:1415,adk:204,fumC:11,gyrB:4,icd:1,mdh:8,purA:8,recA:2"
        #     species = config_sample["properties"].get("Provided_species", None)
        #     {species: script_to_run}

        #     if mlst_string is not None and species is not None:
        #         shell("/srv/data/tools/git.repositories/SSI-scripts/qcscripts/qcfindgene.sh {}".format(species, mlst_string))
        #         shell("touch {output}")
        #     else:
        #         datahandling.log(log_err, "mlst not set or species not set")
        #         shell("touch {output}")

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
        sample = config_sample.get("name", "ERROR") + "__" + component + ".yaml",
    conda:
        "../envs/python_packages.yaml"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_analyzer.py")
