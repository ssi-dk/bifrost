from bifrostlib import datahandling


def extract_contigs_sum_cov(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction("contigs.sum.cov")
    options = sampleComponentObj.get_options()
    yaml = datahandling.load_yaml(file_path)
    for bin_value in options["cov_bin_values"]:
        total_length = 0
        total_depth = 0
        total_contigs = 0
        for contig in yaml["contig_depth"]:
            if yaml["contig_depth"][contig]["coverage"] >= float(bin_value):
                total_length += yaml["contig_depth"][contig]["total_length"]
                total_depth += yaml["contig_depth"][contig]["total_depth"]
                total_contigs += 1
        results[key]["bin_contigs_at_{}x".format(bin_value)] = total_contigs
        results[key]["bin_length_at_{}x".format(bin_value)] = total_length
        results[key]["bin_coverage_at_{}x".format(bin_value)] = float(total_depth / total_length)
        summary["bin_contigs_at_{}x".format(bin_value)] = total_contigs
        summary["bin_length_at_{}x".format(bin_value)] = total_length
        summary["bin_coverage_at_{}x".format(bin_value)] = float(total_depth / total_length)
    return (summary, results)


def extract_contigs_bin_cov(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction("contigs.bin.cov")
    options = sampleComponentObj.get_options()
    results[key] = datahandling.load_yaml(file_path)
    for bin_value in options["cov_bin_values"]:
        summary["raw_length_at_{}x".format(bin_value)] = results[key]["binned_depth"][bin_value - 1]
    return (summary, results)


def extract_bbuk_log(sampleComponentObj):
    import json
    summary, results, file_path, key = sampleComponentObj.start_data_extraction("log/setup__filter_reads_with_bbduk.err.log")
    buffer = datahandling.read_buffer(file_path)
    bbduk_dict = json.loads("{" + buffer.partition("{")[2])
    results[key] = bbduk_dict
    results[key]["filtered_reads_num"] = results[key]["readsIn"] - results[key]["readsRemoved"]
    summary["filtered_reads_num"] = results[key]["filtered_reads_num"]
    return (summary, results)


def extract_quast_report(sampleComponentObj):
    import re
    summary, results, file_path, key = sampleComponentObj.start_data_extraction("quast/report.tsv")
    buffer = datahandling.read_buffer(file_path)
    results[key]["GC"] = float(re.search("GC \(%\)\t([0-9]+[\.]?[0-9]*)", buffer, re.MULTILINE).group(1))
    results[key]["N50"] = int(re.search("N50\t([0-9]+)", buffer, re.MULTILINE).group(1))
    results[key]["N75"] = int(re.search("N75\t([0-9]+)", buffer, re.MULTILINE).group(1))
    results[key]["L50"] = int(re.search("L50\t([0-9]+)", buffer, re.MULTILINE).group(1))
    results[key]["L75"] = int(re.search("L75\t([0-9]+)", buffer, re.MULTILINE).group(1))
    summary["GC"] = results[key]["GC"]
    summary["N50"] = results[key]["N50"]
    return (summary, results)


def extract_contig_variants(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction("contigs.variants")
    yaml = datahandling.load_yaml(file_path)
    results[key] = yaml
    summary["snp_filter_10x_10%"] = yaml["variant_table"][9][9]
    summary["snp_filter_indels"] = yaml["indels"]
    summary["snp_filter_deletions"] = yaml["deletions"]
    return (summary, results)


def extract_contig_stats(sampleComponentObj):
    import re
    summary, results, file_path, key = sampleComponentObj.start_data_extraction("contigs.stats")
    buffer = datahandling.read_buffer(file_path)
    for line in buffer.split("\n"):
        if line.startswith("SN"):
            results[key][re.sub('[^A-Za-z0-9\s]+', '', line.split("\t")[1]).replace(" ", "_").rstrip("_")] = line.split("\t")[2]
    summary["raw_total_sequences"] = results[key]["raw_total_sequences"]
    summary["reads_mapped"] = results[key]["reads_mapped"]
    summary["reads_unmapped"] = results[key]["reads_unmapped"]
    summary["insert_size_average"] = results[key]["insert_size_average"]
    summary["insert_size_standard_deviation"] = results[key]["insert_size_standard_deviation"]
    return (summary, results)


def datadump(sampleComponentObj, log):
    sampleComponentObj.start_data_dump(log=log)
    sampleComponentObj.run_data_dump_on_function(extract_contigs_sum_cov, log=log)
    sampleComponentObj.run_data_dump_on_function(extract_contigs_bin_cov, log=log)
    sampleComponentObj.run_data_dump_on_function(extract_bbuk_log, log=log)
    sampleComponentObj.run_data_dump_on_function(extract_quast_report, log=log)
    sampleComponentObj.run_data_dump_on_function(extract_contig_variants, log=log)
    sampleComponentObj.run_data_dump_on_function(extract_contig_stats, log=log)
    sampleComponentObj.end_data_dump(log=log)


datadump(
    snakemake.params.sampleComponentObj,
    snakemake.log)
