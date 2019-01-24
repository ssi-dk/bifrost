import pkg_resources
from ruamel.yaml import YAML
import os
import re
from bifrostlib import datahandling
import sys

config = datahandling.load_config()


def extract_bbuk_log(file_path, key, data_dict):
    buffer = datahandling.read_buffer(file_path)
    data_dict["results"][key]["input_reads_num"] = int(re.search("Input:\s*([0-9]+)\sreads", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["filtered_reads_num"] = int(re.search("Result:\s*([0-9]+)\sreads", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["input_reads_bases"] = int(re.search("Input:.*?([0-9]+)\sbases", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["filtered_reads_bases"] = int(re.search("Result:.*?([0-9]+)\sbases", buffer, re.MULTILINE).group(1))
    data_dict["summary"]["filtered_reads_num"] = data_dict["results"][key]["filtered_reads_num"]
    return data_dict


def script__datadump_min_read_check(folder, sample):
    folder = str(folder)
    sample = str(sample)
    data_dict = datahandling.load_sample_component(sample)
    data_dict["summary"] = data_dict.get("summary", {})
    data_dict["results"] = data_dict.get("results", {})

    data_dict = datahandling.datadump_template(data_dict, folder, "log/setup__filter_reads_with_bbduk.err.log", extract_bbuk_log)

    datahandling.save_sample_component(data_dict, sample)

    return 0


script__datadump_min_read_check(snakemake.params.folder, snakemake.params.sample)
