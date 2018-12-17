import pkg_resources
import os
import re
from bifrostlib import datahandling
import sys

config = datahandling.load_config()


def extract_debug_report_txt(file_path, key, data_dict):
    buffer = datahandling.read_buffer(file_path)
    buffer = buffer.split("\n")
    values = []
    if len(buffer) > 1:
        header = buffer[0].split("\t")  # remove the header
        for item in buffer:
            values.append(dict(zip(header, item.split("\t"))))
        data_dict["results"][key]["kraken_output"].append(values)
    return data_dict


def script__datadump_resfinder_db(folder, sample, version):
    folder = str(folder)
    sample = str(sample)
    data_dict = datahandling.load_sample_component(sample)
    data_dict["summary"] = data_dict.get("summary", {})
    data_dict["results"] = data_dict.get("results", {})
    data_dict["database"] = config["ariba_resfinder_database"]
    data_dict["version"] = version

    data_dict = datahandling.datadump_template(data_dict, folder, "debug.report.tsv", extract_report_txt)

    datahandling.save_sample_component(data_dict, sample)

    return 0


script__datadump_resfinder_db(snakemake.params.folder, snakemake.params.sample, snakemake.params.version)
