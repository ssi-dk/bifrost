import pkg_resources
from ruamel.yaml import YAML
import os
import re
from bifrostlib import datahandling
import sys

config = datahandling.load_config()


def extract_bracken_txt(file_path, key, data_dict):
    buffer = datahandling.read_buffer(file_path)
    buffer = buffer.split("\n")
    if len(buffer) > 1:
        for i in range(1, len(buffer) - 1 ):  # skip first line as it's header
            data_dict["results"][key]["species_" + str(i) + "_name"] = buffer[i].split("\t")[0]
            data_dict["results"][key]["species_" + str(i) + "_kraken_assigned_reads"] = buffer[i].split("\t")[3]
            data_dict["results"][key]["species_" + str(i) + "_added_reads"] = buffer[i].split("\t")[4]
            data_dict["results"][key]["species_" + str(i) + "_count"] = int(buffer[i].split("\t")[5].strip())
    return data_dict


def extract_kraken_report_bracken_txt(file_path, key, data_dict):
    buffer = datahandling.read_buffer(file_path)
    buffer = buffer.split("\n")
    if len(buffer) > 2:
        data_dict["results"][key]["unclassified_count"] = int(buffer[0].split("\t")[1])
        data_dict["results"][key]["root"] = int(buffer[1].split("\t")[1])
    return data_dict


def extract_kraken_report_txt(file_path, key, data_dict):
    buffer = datahandling.read_buffer(file_path)
    buffer = buffer.split("\n")
    data_dict["results"][key]["kraken_output"] = []
    for item in buffer:
        data_dict["results"][key]["kraken_output"].append([value.strip() for value in item.split("\t")])
    return data_dict


def species_math(data_dict):
    try:
        if "status" not in data_dict["results"]["kraken_report_bracken_txt"] and "status" not in data_dict["results"]["bracken_txt"] and "species_1_count" in data_dict["results"]["bracken_txt"] and "species_2_count" in data_dict["results"]["bracken_txt"]:
            data_dict["summary"]["percent_unclassified"] = data_dict["results"]["kraken_report_bracken_txt"]["unclassified_count"] / data_dict["results"]["kraken_report_bracken_txt"]["unclassified_count"] + data_dict["results"]["kraken_report_bracken_txt"]["root"]
            data_dict["summary"]["percent_classified_species_1"] = data_dict["results"]["bracken_txt"]["species_1_count"] / data_dict["results"]["kraken_report_bracken_txt"]["unclassified_count"] + data_dict["results"]["kraken_report_bracken_txt"]["root"]
            data_dict["summary"]["name_classified_species_1"] = data_dict["results"]["bracken_txt"]["species_1_name"]
            data_dict["summary"]["percent_classified_species_2"] = data_dict["results"]["bracken_txt"]["species_2_count"] / data_dict["results"]["kraken_report_bracken_txt"]["unclassified_count"] + data_dict["results"]["kraken_report_bracken_txt"]["root"]
            data_dict["summary"]["name_classified_species_2"] = data_dict["results"]["bracken_txt"]["species_2_name"]
    except Exception as e:
        print(e)
    return data_dict


def script__datadump_whats_my_species(folder, sample):
    folder = str(folder)
    sample = str(sample)
    data_dict = datahandling.load_sample_component(sample)
    data_dict["summary"] = data_dict.get("summary", {})
    data_dict["results"] = data_dict.get("results", {})

    data_dict = datahandling.datadump_template(data_dict, folder, "bracken.txt", extract_bracken_txt)
    data_dict = datahandling.datadump_template(data_dict, folder, "kraken_report_bracken.txt", extract_kraken_report_bracken_txt)
    data_dict = datahandling.datadump_template(data_dict, folder, "kraken_report.txt", extract_kraken_report_txt)
    data_dict = species_math(data_dict)

    datahandling.save_sample_component(data_dict, sample)

    return 0


script__datadump_whats_my_species(snakemake.params.folder, snakemake.params.sample)
