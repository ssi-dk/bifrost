import os
import sys
import pandas
import Bio.SeqIO
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "../scripts"))
import datahandling

configfile: "../serumqc_config.yaml"
# requires --config R1_reads={read_location},R2_reads={read_location}
sample = config["Sample"]
global_threads = config["threads"]
global_memory_in_GB = config["memory"]

config_sample = datahandling.load_sample(sample)
R1 = config_sample["reads"]["R1"]
R2 = config_sample["reads"]["R2"]

component = "whats_my_species"

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


rule_name = "setup__filter_reads_with_bbduk"
rule setup__filter_reads_with_bbduk:
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
        rules.setup.output.init_file,
        reads = (R1, R2),
    output:
        filtered_reads = temp(rules.setup.params.folder + "/filtered.fastq")
    params:
        adapters = config.get("adapters_fasta", os.path.join(os.path.dirname(workflow.snakefile), "../resources/adapters.fasta"))
    conda:
        "../envs/bbmap.yaml"
    shell:
        "bbduk.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.reads[0]} in2={input.reads[1]} out={output.filtered_reads} ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tbo minbasequality=14 1> {log.out_file} 2> {log.err_file}"


rule_name = "contaminant_check__classify_reads_kraken_minikraken_db"
rule contaminant_check__classify_reads_kraken_minikraken_db:
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
        filtered_reads = rules.setup__filter_reads_with_bbduk.output.filtered_reads,
    output:
        kraken_report = rules.setup.params.folder + "/kraken_report.txt"
    params:
        db = config.get("kraken_database", os.path.join(os.path.dirname(workflow.snakefile), "../resources/kraken_database"))
    conda:
        "../envs/kraken.yaml"
    shell:
        "kraken --threads {threads} -db {params.db} --fastq-input {input.filtered_reads} 2> {log.err_file} | kraken-report -db {params.db} 1> {output.kraken_report}"


rule_name = "contaminant_check__determine_species_bracken_on_minikraken_results"
rule contaminant_check__determine_species_bracken_on_minikraken_results:
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
        kraken_report = rules.contaminant_check__classify_reads_kraken_minikraken_db.output.kraken_report,
    output:
        bracken = rules.setup.params.folder + "/bracken.txt",
        kraken_report_bracken = rules.setup.params.folder + "/kraken_report_bracken.txt"
    params:
        kmer_dist = config.get("kraken_kmer_dist", os.path.join(os.path.dirname(workflow.snakefile), "../resources/kraken_kmer_dist.txt"))
    conda:
        "../envs/bracken.yaml"
    shell:
        """
        est_abundance.py -i {input.kraken_report} -k {params.kmer_dist} -o {output.bracken} 1> {log.out_file} 2> {log.err_file}
        sort -r -t$'\t' -k7 {output.bracken} -o {output.bracken}
        """

rule_name = "species_check__set_species"
rule species_check__set_species:
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
        bracken = rules.contaminant_check__determine_species_bracken_on_minikraken_results.output.bracken,
    output:
        species = rules.setup.params.folder + "/species.txt",
    params:
        sample = sample,
    run:
        try:
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            sample_db = datahandling.load_sample(sample)
            with open(output.species, "w") as species_file:
                df = pandas.read_table(input.bracken)
                # This try-except will avoid a crash when no species is found by bracken.
                try:
                    sample_db["properties"]["detected_species"] = df["name"].iloc[0]
                except IndexError:
                    sample_db["properties"]["detected_species"] = None
                sample_db["properties"]["provided_species"] = sample_db["properties"].get("provided_species",)
                if sample_db["properties"]["provided_species"] is not None:
                    sample_db["properties"]["species"] = sample_db["properties"]["provided_species"]
                else:
                    sample_db["properties"]["species"] = sample_db["properties"]["detected_species"]
            datahandling.save_sample(sample_db, sample)
        except Exception as e:
            datahandling.log(log_err, str(e))

rule_name = "datadump_whats_my_species"
rule datadump_whats_my_species:
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
        species = rules.species_check__set_species.output.species,
    output:
        summary = touch(rules.all.input)
    params:
        sample = config_sample.get("name", "ERROR") + "__" + component + ".yaml",
        folder = rules.setup.params.folder
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_whats_my_species.py")
