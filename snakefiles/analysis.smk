import os
import sys
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "../scripts"))
import datahandling

configfile: os.path.join(os.path.dirname(workflow.snakefile), "../config.yaml")
# requires --config R1_reads={read_location},R2_reads={read_location}
sample = config["Sample"]
global_threads = config["threads"]
global_memory_in_GB = config["memory"]

config_sample = datahandling.load_sample(sample)

R1 = config_sample["sample"]["R1"]
R2 = config_sample["sample"]["R2"]

component = "analysis"

onsuccess:
    print("Workflow complete")
    datahandling.update_sample_component_success(component, sample)

onerror:
    print("Workflow error")
    datahandling.update_sample_component_failure(component, sample)

rule all:
    input:
        component + "/" + component + "_complete"


rule setup:
    output:
        folder = component
    shell:
        "mkdir {output}"


rule_name = "species__checker_and_setter"
rule species__checker:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = component,
        reads = (R1, R2)
    output:
        check_file = touch(component + "/species_set"),
    params:
        species_value = config_sample.get("species", None),
        kraken_db = config["kraken"]["database"]
    run:
        if params.species_value is None:
            shell("kraken --threads {threads} --db {params.kraken_db} {input.reads[0]} {input.reads[1]} | kraken-report --db {params.kraken_db} > kraken_report.txt")
            shell("grep -oPm1 '.*\sS\s[0-9]+\s+(\K.*)' kraken_report.txt > {output.check_file}")
            with open(output.check_file) as species_check:
                species = species_check.readlines()
                if len(species) == 1:
                    config_sample["species"] = species[0]
                    datahandling.save_sample(config_sample, sample)
        else:
            shell("touch {output.checkfile}")


rule_name = "ariba__resfinder"
rule ariba__resfinder:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = component,
        reads = (R1, R2)
    output:
        folder = component + "/ariba_resfinder",
    params:
        database = config["ariba"]["resfinder"]["database"]
    conda:
        "../envs/ariba.yaml"
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
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        contigs = component + "/ariba_resfinder",
    output:
        report = component + "/abricate_on_resfinder_from_ariba.tsv",
    params:
        database = config["abricate"]["resfinder"]["database"],
        db_name = config["abricate"]["resfinder"]["name"],
    conda:
        "../envs/abricate.yaml"
    shell:
        """
        if [[ -e {input.contigs}/assemblies.fa.gz ]] && [[ -n $(gzip -cd {input.contigs}/assemblies.fa.gz | head -c1) ]];
        then abricate --datadir {params.database} --db {params.db_name} {input.contigs}/assemblies.fa.gz > {output.report} 2> {log.err_file};
        else touch {output.report};
        fi;
        """

rule_name = "ariba__plasmidfinder"
rule ariba__plasmidfinder:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = component,
        reads = (R1, R2)
    output:
        folder = component + "/ariba__plasmidfinder",
    params:
        database = config["ariba"]["plasmidfinder"]["database"]
    conda:
        "../envs/ariba.yaml"
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
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = component + "/ariba__plasmidfinder",
    output:
        report = component + "/abricate_on_plasmidfinder_from_ariba.tsv",
    params:
        database = config["abricate"]["plasmidfinder"]["database"],
        db_name = config["abricate"]["plasmidfinder"]["name"],
    conda:
        "../envs/abricate.yaml"
    shell:
        """
        if [[ -e {input.folder}/assemblies.fa.gz ]] && [[ -n $(gzip -cd {input.folder}/assemblies.fa.gz | head -c1) ]];
        then abricate --datadir {params.database} --db {params.db_name} {input.folder}/assemblies.fa.gz > {output.report} 2> {log.err_file};
        else touch {output.report};
        fi;
        """

rule_name = "ariba__mlst"
rule ariba__mlst:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        check_file = component + "/species_set",
        folder = component,
        reads = (R1, R2)
    output:
        folder = component + "/ariba_mlst",
    params:
        species = config_sample.get("species", None)
    conda:
        "../envs/ariba.yaml"
    run:
        if species_db is not None:
            mlst_species_DB = datahandling.get_mlst_species_DB(species)
            if mlst_species is None:
                touch(output.folder)
            else:
                shell("ariba run {} {} {} {} 1> {} 2> {}".format(mlst_species_DB, input.reads[0], input.reads[1], output.folder, log.out_file, log.err_file))
        else:
            touch(output.folder)
        # os.path.join(os.path.dirname(workflow.snakefile), "../scripts/ariba_mlst.py")

rule_name = "datadump_analysis"
rule datadump_analysis:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        component + "/ariba_mlst",
        component + "/abricate_on_resfinder_from_ariba.tsv",
        component + "/abricate_on_plasmidfinder_from_ariba.tsv",
        component,
    output:
        summary = touch(component + "/" + component + "_complete")
    params:
        sample = sample,
    conda:
        "../envs/ariba.yaml"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_analysis.py")
