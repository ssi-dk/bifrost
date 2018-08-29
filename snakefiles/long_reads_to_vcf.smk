# WIP Don't use yet
import os
import sys
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "../scripts"))
import datahandling

# requires --config reads={read_location},reference={ref_location}
sample = config["Sample"]
global_threads = 4
global_memory_in_GB = 20

config_sample = datahandling.load_sample(sample)
reads = config["reads"]
reference = config["reference"]

component = "long_reads_to_vcf"

rule all:
    input:
        component + "/" + component + "_complete"


rule setup:
    output:
        init_file = touch(temp(component + "/" + component + "_initialized")),
    params:
        folder = component


rule_name = "post_assembly__mapping"
rule post_assembly__mapping:
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
        folderinit = rules.setup.output.init_file,
        reads = reads,
        reference = reference,
    output:
        mapped = temp(rules.setup.params.folder + "/contigs.sam")
    shell:
        "minimap2 -t {threads} --MD -ax map-ont {input.reference} {input.reads} 1> {output.mapped} 2> {log.err_file}"


rule_name = "sam_to_sorted_bam"
rule sam_to_sorted_bam:
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
        mapped = rules.post_assembly__mapping.output.mapped,
    output:
        unsorted_bam = temp(rules.setup.params.folder + "/contigs.unsorted.bam"),
        bam = rules.setup.params.folder + "/contigs.bam",
    conda:
        "../envs/minimap2.yaml"
    shell:
        # Be aware of samtools version
        """
        samtools view -S -b {input.mapped} > {output.unsorted_bam}
        samtools sort -@{threads} -o {output.bam} {output.unsorted_bam}
        """

rule_name = "post_assembly__pileup"
rule post_assembly__pileup:
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
        mapped = rules.post_assembly__mapping.output.mapped
    output:
        coverage = temp(rules.setup.params.folder + "/contigs.cov"),
        pileup = rules.setup.params.folder + "/contigs.pileup"
    conda:
        "../envs/bbmap.yaml"
    shell:
        "pileup.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.mapped} basecov={output.coverage} out={output.pileup} 1> {log.out_file} 2> {log.err_file}"


rule_name = "post_assembly__call_variants"
rule post_assembly__call_variants:
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
        contigs = reference,
        mapped = rules.sam_to_sorted_bam.output.bam
    output:
        variants = rules.setup.params.folder + "/contigs.vcf",
        final = touch(rules.all.input)
    conda:
        "../envs/bbmap.yaml"
    shell:
        "callvariants.sh in={input.mapped} vcf={output.variants} ref={input.contigs} ploidy=1 clearfilters 1> {log.out_file} 2> {log.err_file}"
