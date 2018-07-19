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

# my understanding is all helps specify final output
component = "assembly"

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
        folder = directory(component)
    shell:
        "mkdir {output}"

rule setup__filter_reads_with_bbduk:
    message:
        "Running step: {rule}"
    input:
        directory = "assembly",
        reads = (R1, R2)
    output:
        filtered_reads = temp("assembly/filtered.fastq")
    params:
        adapters = os.path.join(os.path.dirname(workflow.snakefile), "../resources/adapters.fasta")
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/bbmap.yaml"
    log:
        "assembly/log/setup__filter_reads_with_bbduk.log"
    benchmark:
        "assembly/benchmarks/setup__filter_reads_with_bbduk.benchmark"
    shell:
        "bbduk.sh in={input.reads[0]} in2={input.reads[1]} out={output.filtered_reads} ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tbo minbasequality=14 &> {log}"


# rule assembly_check__combine_reads_with_bbmerge:
#     message:
#         "Running step: {rule}"
#     input:
#         filtered_reads = "assembly/filtered.fastq"
#     output:
#         merged_reads = "assembly/merged.fastq",
#         unmerged_reads = "assembly/unmerged.fastq"
#     threads:
#         global_threads
#     resources:
#         memory_in_GB = global_memory_in_GB
#     log:
#         "assembly/log/assembly_check__combine_reads_with_bbmerge.log"
#     benchmark:
#         "assembly/benchmarks/assembly_check__combine_reads_with_bbmerge.benchmark"
#     run:
#         if os.path.isfile("qcquickie/merged.fastq") and os.path.isfile("qcquickie/unmerged.fastq"):
#             shell("ln -s {} {}".format(os.path.realpath("qcquickie/merged.fastq"), output.merged_reads))
#             shell("ln -s {} {}".format(os.path.realpath("qcquickie/unmerged.fastq"), output.unmerged_reads))
#         else:
#             shell("bbmerge.sh in={input.filtered_reads} out={output.merged_reads} outu={output.unmerged_reads} &> {log}")


rule assembly__spades:
    message:
        "Running step: {rule}"
    input:
        filtered_reads = "assembly/filtered.fastq",
    params:
        assembler = "assembler/assembler_SPAdes"
    output:
        spades_folder = temp("spades"),
        contigs = "assembly/temp.fasta",
        assembly_with = touch("assembly/assembly_with_SPAdes"),
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/spades.yaml"
    log:
        "assembly/log/assembly__spades.log"
    benchmark:
        "assembly/benchmarks/assembly__spades.benchmark"
    shell:
        """
        spades.py -k 21,33,55,77 --12 {input.filtered_reads} -o {output.spades_folder} --careful &> {log}
        mv {spades_folder}/contigs.fasta {output.contigs}
        """

rule assembly__skesa:
    message:
        "Running step: {rule}"
    input:
        filtered_reads = "assembly/filtered.fastq",
    output:
        contigs = "assembly/temp.fasta",
        assembly_with = touch("assembly/assembly_with_skesa")
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/skesa.yaml"
    log:
        "assembly/log/assembly__skesa.log"
    benchmark:
        "assembly/benchmarks/assembly__skesa.benchmark"
    shell:
        "skesa --cores {threads} --memory {resources.memory_in_GB} --use_paired_ends --fastq {input.filtered_reads} --contigs_out {output.contigs} &> {log}"


rule assembly__selection:
    input:
        assembly_with = "assembly/assembly_with_" + config["assembly_with"]
    params:
        "assembly/temp.fasta"
    output:
        "assembly/contigs.fasta"
    shell:
        "mv {params} {output}"


rule assembly_check__quast_on_contigs:
    message:
        "Running step: {rule}"
    input:
        contigs = "assembly/contigs.fasta"
    output:
        quast = directory("assembly/quast")
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/quast.yaml"
    log:
        "assembly/log/assembly_check__quast_on_tadpole_contigs.log"
    benchmark:
        "assembly/benchmarks/assembly_check__quast_on_tadpole_contigs.benchmark"
    shell:
        "quast.py --threads {threads} {input.contigs} -o {output.quast} &> {log}"


rule assembly_check__sketch_on_contigs:
    message:
        "Running step: {rule}"
    input:
        contigs = "assembly/contigs.fasta"
    output:
        sketch = "assembly/contigs.sketch"
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/bbmap.yaml"
    log:
        "assembly/log/assembly_check__sketch_on_contigs.log"
    benchmark:
        "assembly/benchmarks/assembly_check__sketch_on_contigs.benchmark"
    shell:
        "sketch.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.contigs} out={output.sketch} &> {log}"


rule post_assembly__stats:
    message:
        "Running step: {rule}"
    input:
        contigs = "assembly/contigs.fasta"
    output:
        stats = touch("assembly/post_assermbly__stats")
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/bbmap.yaml"
    log:
        "assembly/log/post_assembly__stats.log"
    benchmark:
        "assembly/benchmarks/post_assembly__stats.benchmark"
    shell:
        "stats.sh {input.contigs} &> {log}"


# rule post_assembly__mapping:
#     message:
#         "Running step: {rule}"
#     input:
#         contigs = "assembly/contigs.fasta",
#         filtered_reads = "assembly/filtered.fastq",
#     output:
#         mapped = "assembly/contigs.sam"
#     threads:
#         global_threads
#     resources:
#         memory_in_GB = global_memory_in_GB
#     conda:
#         "../envs/bbmap.yaml"
#     group:
#         "assembly"
#     log:
#         "assembly/log/post_assembly__mapping.log"
#     benchmark:
#         "assembly/benchmarks/post_assembly__mapping.benchmark"
#     shell:
#         "bbmap.sh threads={threads} -Xmx{resources.memory_in_GB}G ref={input.contigs} in={input.filtered_reads} out={output.mapped} &> {log}"


rule post_assembly__mapping:
    message:
        "Running step: {rule}"
    input:
        contigs = "assembly/contigs.fasta",
        filtered_reads = "assembly/filtered.fastq",
    output:
        mapped = temp("assembly/contigs.sam")
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/minimap2.yaml"
    log:
        "assembly/log/post_assembly__mapping.log"
    benchmark:
        "assembly/benchmarks/post_assembly__mapping.benchmark"
    shell:
        "minimap2 -t {threads} --MD -ax sr {input.contigs} {input.filtered_reads} 1> {output.mapped} 2> {log}"


rule post_assembly__samtools_stats:
    message:
        "Running step: {rule}"
    input:
        mapped = "assembly/contigs.sam"
    output:
        stats = "assembly/contigs.stats",
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/samtools.yaml"
    log:
        "assembly/log/post_assembly__samtools_stats.log"
    benchmark:
        "assembly/benchmarks/post_assembly__samtools_stats.benchmark"
    shell:
        "samtools stats -@ {threads} {input.mapped} 1> {output.stats} 2> {log}"


rule post_assembly__pileup:
    message:
        "Running step: {rule}"
    input:
        mapped = "assembly/contigs.sam"
    output:
        coverage = temp("assembly/contigs.cov"),
        pileup = "assembly/contigs.pileup"
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/bbmap.yaml"
    log:
        "assembly/log/post_assembly__pileup.log"
    benchmark:
        "assembly/benchmarks/post_assembly__pileup.benchmark"
    shell:
        "pileup.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.mapped} basecov={output.coverage} out={output.pileup} &> {log}"


rule summarize__depth:
    message:
        "Running step: {rule}"
    input:
        coverage = "assembly/contigs.cov"
    output:
        contig_depth_yaml = "assembly/contigs.sum.cov",
        binned_depth_yaml = "assembly/contigs.bin.cov"
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/python_packages.yaml"
    log:
        "assembly/log/summarize__bin_coverage.log"
    benchmark:
        "assembly/benchmarks/summarize__bin_coverage.benchmark"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/summarize_depth.py")


rule post_assembly__call_variants:
    message:
        "Running step: {rule}"
    input:
        contigs = "assembly/contigs.fasta",
        mapped = "assembly/contigs.sam",
    output:
        variants = temp("assembly/contigs.vcf"),
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/bbmap.yaml"
    log:
        "assembly/log/post_assembly__call_variants.log"
    benchmark:
        "assembly/benchmarks/post_assembly__call_variants.benchmark"
    shell:
        "callvariants.sh in={input.mapped} vcf={output.variants} ref={input.contigs} ploidy=1 clearfilters &> {log}"


rule summarize__variants:
    message:
        "Running step: {rule}"
    input:
        variants = "assembly/contigs.vcf",
    output:
        variants_yaml = "assembly/contigs.variants",
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/python_packages.yaml"
    log:
        "assembly/log/summarize__variants.log"
    benchmark:
        "assembly/benchmarks/summarize__variants.benchmark"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/summarize_variants.py")


rule post_assembly__annotate:
    message:
        "Running step: {rule}"
    input:
        contigs = "assembly/contigs.fasta"
    output:
        gff = "assembly/contigs.gff",
    params:
        prokka = temp("assembly/prokka")
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "../envs/prokka.yaml"
    log:
        "assembly/log/post_assembly__annotate.log"
    benchmark:
        "assembly/benchmarks/post_assembly__annotate.benchmark"
    shell:
        """ 
        prokka --cpus {threads} --centre XXX --compliant --outdir {params.prokka} {input.contigs} &> {log} 
        mv {params.prokka}/*.gff {output.gff} 
        """ 

rule datadump_assembly:
    message:
        "Running step: {rule}"
    input:
        "assembly/contigs.gff",
        "assembly/contigs.bin.cov",
        "assembly/contigs.sum.cov",
        "assembly/contigs.variants",
        "assembly/quast",
        "assembly/contigs.stats",
        "assembly/contigs.sketch",
        folder = "assembly",
    output:
        summary = touch(rules.all.input)
    params:
        sample = sample,
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        "assembly/log/datadump_assembly.log"
    benchmark:
        "assembly/benchmarks/datadump_assembly.benchmark"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_assembly.py")
