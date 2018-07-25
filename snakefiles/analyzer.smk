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

R1 = config_sample["reads"]["R1"]
R2 = config_sample["reads"]["R2"]

component = "analyzer"


onsuccess:
    print("Workflow complete")
    datahandling.update_sample_component_success(config_sample.get("name","ERROR") + "__" + component + ".yaml")


onerror:
    print("Workflow error")
    datahandling.update_sample_component_failure(config_sample.get("name","ERROR") + "__" + component + ".yaml")


rule all:
    input:
        rules.setup.output.folder + "/" + component + "_complete"


rule setup:
    output:
        folder = directory(component)
    shell:
        "mkdir {output}"


rule_name = "species_checker_and_setter"
rule species_checker_and_setter:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = rules.setup.output.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.output.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.output.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = rules.setup.output.folder,
        reads = (R1, R2)
    output:
        check_file = touch(rules.setup.output.folder + "/species_set"),
    params:
        kraken_db = config["kraken"]["database"]
    run:
        sample_db = datahandling.load_sample(sample)
        if sample_db["properties"].get("species",) is None:
            shell("kraken --threads {threads} --db {params.kraken_db} {input.reads[0]} {input.reads[1]} | kraken-report --db {params.kraken_db} > kraken_report.txt")
            shell("grep -oPm1 '.*\sS\s[0-9]+\s+(\K.*)' kraken_report.txt > {output.check_file}")
            with open(output.check_file) as species_check:
                species = species_check.readlines()
                if len(species) == 1:
                    sample_db["properties"]["species"] = species[0].strip()
                    datahandling.save_sample(sample_db, sample)
        else:
            shell("touch {output.check_file}")


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
        out_file = rules.setup.output.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.output.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.output.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = rules.setup.output.folder,
        reads = (R1, R2)
    output:
        folder = directory(rules.setup.output.folder + "/ariba_resfinder")
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
        out_file = rules.setup.output.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.output.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.output.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        contigs = rules.ariba_resfinder.output.folder
    output:
        report = rules.setup.output.folder + "/abricate_on_resfinder_from_ariba.tsv",
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


rule_name = "ariba_plasmidfinder"
rule ariba_plasmidfinder:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = rules.setup.output.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.output.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.output.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = rules.setup.output.folder,
        reads = (R1, R2)
    output:
        folder = directory(rules.setup.output.folder + "/ariba_plasmidfinder")
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
        out_file = rules.setup.output.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.output.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.output.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = rules.ariba_plasmidfinder.output
    output:
        report = rules.setup.output.folder + "/abricate_on_plasmidfinder_from_ariba.tsv",
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
        out_file = rules.setup.output.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.output.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.output.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        check_file = rules.species_checker_and_setter.output,
        folder = rules.setup.output.folder,
        reads = (R1, R2)
    output:
        folder = directory(rules.setup.output.folder + "/ariba_mlst")
    params:
        sample = sample
    conda:
        "../envs/ariba.yaml"
    run:
        mlst_species_DB = datahandling.get_mlst_species_DB(sample)
        if mlst_species_DB is None:
            touch(output.folder)
        else:
            shell("ariba run {} {} {} {} 1> {} 2> {}".format(mlst_species_DB, input.reads[0], input.reads[1], output.folder, log.out_file, log.err_file))


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
        out_file = rules.setup.output.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.output.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.output.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        rules.ariba_mlst.output.folder,
        rules.abricate_on_ariba_resfinder.output.report,
        rules.abricate_on_ariba_plasmidfinder.output.report,
    output:
        summary = touch(rules.all.input)
    params:
        sample = config_sample.get("name","ERROR") + "__" + component + ".yaml",
    conda:
        "../envs/ariba.yaml"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_analyzer.py")
