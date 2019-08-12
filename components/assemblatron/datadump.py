import pkg_resources
from ruamel.yaml import YAML
import os
import re
from bifrostlib import datahandling
import sys

config = datahandling.load_config()

global GLOBAL_BIN_VALUES
GLOBAL_BIN_VALUES = [1, 10, 25]


def extract_contigs_sum_cov(file_path, key, data_dict):
    yaml = datahandling.load_yaml(file_path)
    data_dict["results"][key] = yaml
    for bin_value in GLOBAL_BIN_VALUES:
        total_length = 0
        total_depth = 0
        total_contigs = 0
        for contig in yaml["contig_depth"]:
            if yaml["contig_depth"][contig]["coverage"] >= float(bin_value):
                total_length += yaml["contig_depth"][contig]["total_length"]
                total_depth += yaml["contig_depth"][contig]["total_depth"]
                total_contigs += 1
        data_dict["summary"]["bin_contigs_at_{}x".format(bin_value)] = total_contigs
        data_dict["summary"]["bin_length_at_{}x".format(bin_value)] = total_length
        data_dict["summary"]["bin_coverage_at_{}x".format(bin_value)] = float(total_depth / total_length)
    return data_dict


def extract_contigs_bin_cov(file_path, key, data_dict):
    yaml = datahandling.load_yaml(file_path)
    data_dict["results"][key] = yaml
    for bin_value in GLOBAL_BIN_VALUES:
        data_dict["summary"]["raw_length_at_{}x".format(bin_value)] = yaml["binned_depth"][bin_value - 1]
    return data_dict


def extract_bbuk_log(file_path, key, data_dict):
    buffer = datahandling.read_buffer(file_path)
    data_dict["results"][key]["input_reads_num"] = int(re.search("Input:\s*([0-9]+)\sreads", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["filtered_reads_num"] = int(re.search("Result:\s*([0-9]+)\sreads", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["input_reads_bases"] = int(re.search("Input:.*?([0-9]+)\sbases", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["filtered_reads_bases"] = int(re.search("Result:.*?([0-9]+)\sbases", buffer, re.MULTILINE).group(1))
    data_dict["summary"]["filtered_reads_num"] = data_dict["results"][key]["filtered_reads_num"]
    return data_dict


def extract_quast_report(file_path, key, data_dict):
    buffer = datahandling.read_buffer(file_path)
    data_dict["results"][key]["GC"] = float(re.search("GC \(%\)\t([0-9]+[\.]?[0-9]*)", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["N50"] = int(re.search("N50\t([0-9]+)", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["N75"] = int(re.search("N75\t([0-9]+)", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["L50"] = int(re.search("L50\t([0-9]+)", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["L75"] = int(re.search("L75\t([0-9]+)", buffer, re.MULTILINE).group(1))
    data_dict["summary"]["GC"] = data_dict["results"][key]["GC"]
    data_dict["summary"]["N50"] = data_dict["results"][key]["N50"]
    return data_dict


def extract_contig_variants(file_path, key, data_dict):
    yaml = datahandling.load_yaml(file_path)
    data_dict["results"][key] = yaml
    data_dict["summary"]["snp_filter_10x_10%"] = yaml["variant_table"][9][9]
    data_dict["summary"]["snp_filter_indels"] = yaml["indels"]
    data_dict["summary"]["snp_filter_deletions"] = yaml["deletions"]
    return data_dict


def extract_contig_stats(file_path, key, data_dict):
    buffer = datahandling.read_buffer(file_path)
    for line in buffer.split("\n"):
        if line.startswith("SN"):
            data_dict["results"][key][re.sub('[^A-Za-z0-9\s]+', '', line.split("\t")[1]).replace(" ", "_").rstrip("_")] = line.split("\t")[2]
    data_dict["summary"]["raw_total_sequences"] = data_dict["results"][key]["raw_total_sequences"]
    data_dict["summary"]["reads_mapped"] = data_dict["results"][key]["reads_mapped"]
    data_dict["summary"]["reads_unmapped"] = data_dict["results"][key]["reads_unmapped"]
    data_dict["summary"]["insert_size_average"] = data_dict["results"][key]["insert_size_average"]
    data_dict["summary"]["insert_size_standard_deviation"] = data_dict["results"][key]["insert_size_standard_deviation"]
    return data_dict


def extract_contig_sketch(file_path, key, data_dict):
    buffer = datahandling.read_buffer(file_path)
    data_dict["results"][key] = buffer.split("\n")
    return data_dict


def script__datadump(folder, sample_file, component_file, sample_component_file, log):
    try:
        log_out = str(log.out_file)
        log_err = str(log.err_file)
        db_sample = datahandling.load_sample(sample_file)
        db_component = datahandling.load_component(component_file)
        db_sample_component = datahandling.load_sample_component(sample_component_file)
        this_function_name = sys._getframe().f_code.co_name

        datahandling.log(log_out, "Started {}\n".format(this_function_name))

        # Initialization of values, summary and reporter are also saved into the sample
        db_sample_component["summary"] = {"component": {"id": db_component["_id"], "date": datetime.datetime.utcnow()}}
        db_sample_component["results"] = {}

        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "contigs.sum.cov", extract_contigs_sum_cov)
        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "contigs.bin.cov", extract_contigs_bin_cov)
        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "log/setup__filter_reads_with_bbduk.err.log", extract_bbuk_log)
        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "quast/report.tsv", extract_quast_report)
        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "contigs.variants", extract_contig_variants)
        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "contigs.stats", extract_contig_stats)
        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "contigs.sketch", extract_contig_sketch)

        datahandling.save_sample_component(db_sample_component, sample_component_file)

    except Exception:
        datahandling.log(log_out, "Exception in {}\n".format(this_function_name))
        datahandling.log(log_err, str(traceback.format_exc()))

    finally:
        datahandling.log(log_out, "Done {}\n".format(this_function_name))
        return 0


script__datadump(
    snakemake.params.folder,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.sample_component_file,
    snakemake.log)
