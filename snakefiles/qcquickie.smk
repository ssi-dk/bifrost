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

component = "qcquickie"

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


rule_name = "fastqc_on_reads"
rule fastqc_on_reads:
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
        directory = rules.setup.output.init_file,
        reads = (R1, R2)
    output:
        folder = directory(rules.setup.params.folder + "/fastqc"),
        fastqc_summary = rules.setup.params.folder + "/fastqc_data.txt"
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mkdir {output.folder}
        fastqc --extract -o {output.folder} -t {threads} {input.reads[0]} {input.reads[1]} 1> {log.out_file} 2> {log.err_file}
        cat {output.folder}/*/fastqc_data.txt > {output.fastqc_summary}
        """


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
        directory = rules.setup.output.init_file,
        reads = (R1, R2)
    output:
        filtered_reads = temp(rules.setup.params.folder + "/filtered.fastq")
    params:
        adapters = config.get("adapters_fasta", os.path.join(os.path.dirname(workflow.snakefile), "../resources/adapters.fasta"))
    conda:
        "../envs/bbmap.yaml"
    shell:
        "bbduk.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.reads[0]} in2={input.reads[1]} out={output.filtered_reads} ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tbo minbasequality=14 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly_check__combine_reads_with_bbmerge"
rule assembly_check__combine_reads_with_bbmerge:
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
        merged_reads = temp(rules.setup.params.folder + "/merged.fastq"),
        unmerged_reads = temp(rules.setup.params.folder + "/unmerged.fastq")
    conda:
        "../envs/bbmap.yaml"
    shell:
        "bbmerge.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.filtered_reads} out={output.merged_reads} outu={output.unmerged_reads} 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly_check__quick_assembly_with_tadpole"
rule assembly_check__quick_assembly_with_tadpole:
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
        merged_reads = rules.assembly_check__combine_reads_with_bbmerge.output.merged_reads,
        unmerged_reads = rules.assembly_check__combine_reads_with_bbmerge.output.unmerged_reads
    output:
        contigs = temp(rules.setup.params.folder + "/raw_contigs.fasta")
    conda:
        "../envs/bbmap.yaml"
    shell:
        "tadpole.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.merged_reads},{input.unmerged_reads} out={output.contigs} 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly_check__rename_contigs"
rule assembly_check__rename_contigs:
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
        contigs = rules.assembly_check__quick_assembly_with_tadpole.output.contigs,
    output:
        contigs = rules.setup.params.folder + "/contigs.fasta"
    run:
        try:
            input_file = str(input.contigs)
            output_file = str(output.contigs)
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            with open(input_file, "r") as fasta_input:
                records = list(Bio.SeqIO.parse(fasta_input, "fasta"))
            for record in records:
                record.id = record.id.split(",")[0]
                record.description = record.id
            with open(output_file, "w") as output_handle:
                Bio.SeqIO.write(records, output_handle, "fasta")
            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(e))


rule_name = "assembly_check__quast_on_contigs"
rule assembly_check__quast_on_contigs:
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
    message:
        "Running step: {rule}"
    input:
        contigs = rules.assembly_check__rename_contigs.output.contigs
    output:
        quast = directory(rules.setup.params.folder + "/quast")
    conda:
        "../envs/quast.yaml"
    shell:
        "quast.py --threads {threads} {input.contigs} -o {output.quast} 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly_check__sketch_on_contigs"
rule assembly_check__sketch_on_contigs:
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
        contigs = rules.assembly_check__rename_contigs.output.contigs
    output:
        sketch = rules.setup.params.folder + "/contigs.sketch"
    conda:
        "../envs/bbmap.yaml"
    shell:
        "sketch.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.contigs} out={output.sketch} 1> {log.out_file} 2> {log.err_file}"


rule_name = "post_assembly__stats"
rule post_assembly__stats:
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
    message:
        "Running step: {rule}"
    input:
        contigs = rules.assembly_check__rename_contigs.output.contigs
    output:
        stats = touch(rules.setup.params.folder + "/post_assermbly__stats")
    conda:
        "../envs/bbmap.yaml"
    shell:
        "stats.sh {input.contigs} 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly_check__map_reads_to_assembly_with_bbmap"
rule assembly_check__map_reads_to_assembly_with_bbmap:
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
        contigs = rules.assembly_check__rename_contigs.output.contigs,
        filtered = rules.setup__filter_reads_with_bbduk.output.filtered_reads
    output:
        mapped = temp(rules.setup.params.folder + "/contigs.sam")
    conda:
        "../envs/bbmap.yaml"
    shell:
        "bbmap.sh threads={threads} -Xmx{resources.memory_in_GB}G ref={input.contigs} in={input.filtered} out={output.mapped} ambig=random 1> {log.out_file} 2> {log.err_file}"


rule_name = "post_assembly__samtools_stats"
rule post_assembly__samtools_stats:
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
        mapped = rules.assembly_check__map_reads_to_assembly_with_bbmap.output.mapped
    output:
        stats = rules.setup.params.folder + "/contigs.stats",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats -@ {threads} {input.mapped} 1> {output.stats} 2> {log.err_file}"


rule_name = "assembly_check__pileup_on_mapped_reads"
rule assembly_check__pileup_on_mapped_reads:
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
        mapped = rules.assembly_check__map_reads_to_assembly_with_bbmap.output.mapped
    output:
        coverage = temp(rules.setup.params.folder + "/contigs.cov"),
        pileup = rules.setup.params.folder + "/contigs.pileup"
    conda:
        "../envs/bbmap.yaml"
    shell:
        "pileup.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.mapped} basecov={output.coverage} out={output.pileup} 1> {log.out_file} 2> {log.err_file}"


rule_name = "summarize__depth"
rule summarize__depth:
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
        coverage = rules.assembly_check__pileup_on_mapped_reads.output.coverage
    output:
        contig_depth_yaml = rules.setup.params.folder + "/contigs.sum.cov",
        binned_depth_yaml = rules.setup.params.folder + "/contigs.bin.cov"
    conda:
        "../envs/python_packages.yaml"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/summarize_depth.py")


rule_name = "assembly_check__call_variants"
rule assembly_check__call_variants:
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
        contigs = rules.assembly_check__rename_contigs.output.contigs,
        mapped = rules.assembly_check__map_reads_to_assembly_with_bbmap.output.mapped,
    output:
        variants = temp(rules.setup.params.folder + "/contigs.vcf")
    conda:
        "../envs/bbmap.yaml"
    shell:
        "callvariants.sh in={input.mapped} vcf={output.variants} ref={input.contigs} ploidy=1 clearfilters 1> {log.out_file} 2> {log.err_file}"


rule_name = "summarize__variants"
rule summarize__variants:
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
        variants = rules.assembly_check__call_variants.output.variants,
    output:
        variants_yaml = rules.setup.params.folder + "/contigs.variants",
    conda:
        "../envs/python_packages.yaml"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/summarize_variants.py")


rule_name = "datadump_qcquickie"
rule datadump_qcquickie:
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
        rules.fastqc_on_reads.output.fastqc_summary,
        rules.assembly_check__quast_on_contigs.output.quast,
        rules.assembly_check__sketch_on_contigs.output.sketch,
        rules.post_assermbly__stats.output.stats,
        rules.post_assembly__samtools_stats.output.stats,
        rules.summarize__depth.output.contig_depth_yaml,
        rules.summarize__variants.output.variants_yaml,
    output:
        summary = touch(rules.all.input)
    params:
        sample = config_sample.get("name", "ERROR") + "__" + component + ".yaml",
        folder = rules.setup.params.folder,
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_qcquickie.py")
