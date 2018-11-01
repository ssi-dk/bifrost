import os
import sys
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "../scripts"))
import datahandling

configfile: "../run_config.yaml"
# requires --config R1_reads={read_location},R2_reads={read_location}
sample = config["Sample"]
global_threads = config["threads"]
global_memory_in_GB = config["memory"]

config_sample = datahandling.load_sample(sample)
R1 = config_sample["reads"]["R1"]
R2 = config_sample["reads"]["R2"]

component = "assemblatron"


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
        folderinit = rules.setup.output.init_file,
        reads = (R1, R2)
    output:
        filtered_reads = temp(rules.setup.params.folder + "/filtered.fastq")
    params:
        adapters = os.path.join(os.path.dirname(workflow.snakefile), config["adapters_fasta"])
    conda:
        "../envs/bbmap.yaml"
    shell:
        "bbduk.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.reads[0]} in2={input.reads[1]} out={output.filtered_reads} ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tbo minbasequality=14 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly__spades"
rule assembly__spades:
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
        spades_folder = temp(directory("spades")),
        contigs = rules.setup.params.folder + "/temp.fasta",
        assembly_with = touch(rules.setup.params.folder + "/assembly_with_SPAdes"),
    conda:
        "../envs/spades.yaml"
    shell:
        """
        spades.py -k 21,33,55,77 --12 {input.filtered_reads} -o {output.spades_folder} --careful 1> {log.out_file} 2> {log.err_file}
        mv {spades_folder}/contigs.fasta {output.contigs}
        """


rule_name = "assembly__skesa"
rule assembly__skesa:
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
        contigs = rules.setup.params.folder + "/temp.fasta",
        assembly_with = touch(rules.setup.params.folder + "/assembly_with_skesa"),
    conda:
        "../envs/skesa.yaml"
    shell:
        "skesa --cores {threads} --memory {resources.memory_in_GB} --use_paired_ends --fastq {input.filtered_reads} --contigs_out {output.contigs} 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly__selection"
rule assembly__selection:
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
        assembly_with = rules.setup.params.folder + "/assembly_with_" + config["assembly_with"]
    params:
        rules.setup.params.folder + "/temp.fasta"
    output:
        rules.setup.params.folder + "/contigs.fasta"
    shell:
        "mv {params} {output}"

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
    input:
        contigs = rules.assembly__selection.output
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
        contigs = rules.assembly__selection.output
    output:
        sketch = rules.setup.params.folder + "/contigs.sketch"
    conda:
        "../envs/bbmap.yaml"
    shell:
        "sendsketch.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.contigs} outsketch={output.sketch} 1> {log.out_file} 2> {log.err_file}"


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
        contigs = rules.assembly__selection.output
    output:
        stats = touch(rules.setup.params.folder + "/post_assermbly__stats")
    conda:
        "../envs/bbmap.yaml"
    shell:
        "stats.sh threads={threads} -Xmx{resources.memory_in_GB}G {input.contigs} 1> {log.out_file} 2> {log.err_file}"


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
        contigs = rules.assembly__selection.output,
        filtered_reads = rules.setup__filter_reads_with_bbduk.output.filtered_reads
    output:
        mapped = temp(rules.setup.params.folder + "/contigs.sam")
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -t {threads} --MD -ax sr {input.contigs} {input.filtered_reads} 1> {output.mapped} 2> {log.err_file}"


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
        mapped = rules.post_assembly__mapping.output.mapped
    output:
        stats = rules.setup.params.folder + "/contigs.stats",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats -@ {threads} {input.mapped} 1> {output.stats} 2> {log.err_file}"


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
        coverage = rules.post_assembly__pileup.output.coverage
    output:
        contig_depth_yaml = rules.setup.params.folder + "/contigs.sum.cov",
        binned_depth_yaml = rules.setup.params.folder + "/contigs.bin.cov"
    conda:
        "../envs/python_packages.yaml"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/summarize_depth.py")


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
        contigs = rules.assembly__selection.output,
        mapped = rules.post_assembly__mapping.output.mapped
    output:
        variants = temp(rules.setup.params.folder + "/contigs.vcf"),
    conda:
        "../envs/bbmap.yaml"
    shell:
        "callvariants.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.mapped} vcf={output.variants} ref={input.contigs} ploidy=1 clearfilters 1> {log.out_file} 2> {log.err_file}"


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
        variants = rules.post_assembly__call_variants.output.variants
    output:
        variants_yaml = rules.setup.params.folder + "/contigs.variants",
    conda:
        "../envs/python_packages.yaml"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/summarize_variants.py")


rule_name = "post_assembly__annotate"
rule post_assembly__annotate:
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
        contigs = rules.assembly__selection.output,
    output:
        gff = rules.setup.params.folder + "/contigs.gff",
    params:
        prokka = temp(directory(rules.setup.params.folder + "/prokka"))
    conda:
        "../envs/prokka.yaml"
    shell:
        """ 
        prokka --cpus {threads} --centre XXX --compliant --outdir {params.prokka} {input.contigs} 1> {log.out_file} 2> {log.err_file};
        mv {params.prokka}/*.gff {output.gff};
        """ 


rule_name = "datadump_assemblatron"
rule datadump_assemblatron:
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
        rules.post_assembly__annotate.output.gff,
        rules.summarize__depth.output.contig_depth_yaml,
        rules.summarize__depth.output.binned_depth_yaml,
        rules.summarize__variants.output.variants_yaml,
        rules.assembly_check__quast_on_contigs.output.quast,
        rules.post_assembly__stats.output.stats,
        rules.post_assembly__samtools_stats.output.stats,
        rules.assembly_check__sketch_on_contigs.output.sketch,
    output:
        summary = touch(rules.all.input)
    params:
        sample = config_sample.get("name", "ERROR") + "__" + component + ".yaml",
        folder = rules.setup.params.folder,
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "../scripts/datadump_assemblatron.py")
