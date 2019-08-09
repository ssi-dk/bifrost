import pkg_resources
import datetime
import os
import re
import sys
import traceback
from bifrostlib import datahandling

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


def species_math(file_path, key, data_dict):
    if "status" not in data_dict["results"]["kraken_report_bracken_txt"] and "status" not in data_dict["results"]["bracken_txt"] and "species_1_count" in data_dict["results"]["bracken_txt"] and "species_2_count" in data_dict["results"]["bracken_txt"]:
        data_dict["summary"]["percent_unclassified"] = data_dict["results"]["kraken_report_bracken_txt"]["unclassified_count"] / (data_dict["results"]["kraken_report_bracken_txt"]["unclassified_count"] + data_dict["results"]["kraken_report_bracken_txt"]["root"])
        data_dict["summary"]["percent_classified_species_1"] = data_dict["results"]["bracken_txt"]["species_1_count"] / (data_dict["results"]["kraken_report_bracken_txt"]["unclassified_count"] + data_dict["results"]["kraken_report_bracken_txt"]["root"])
        data_dict["summary"]["name_classified_species_1"] = data_dict["results"]["bracken_txt"]["species_1_name"]
        data_dict["summary"]["percent_classified_species_2"] = data_dict["results"]["bracken_txt"]["species_2_count"] / (data_dict["results"]["kraken_report_bracken_txt"]["unclassified_count"] + data_dict["results"]["kraken_report_bracken_txt"]["root"])
        data_dict["summary"]["name_classified_species_2"] = data_dict["results"]["bracken_txt"]["species_2_name"]
    return data_dict


def set_sample_species(file_path, key, data_dict):
    if data_dict["properties"]["provided_species"] is not None:
        data_dict["properties"]["species"] = data_dict["properties"]["provided_species"]
    else:
        data_dict["properties"]["species"] = data_dict["properties"]["detected_species"]
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

        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "bracken.txt", extract_bracken_txt)
        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "kraken_report_bracken.txt", extract_kraken_report_bracken_txt)
        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "kraken_report.txt", extract_kraken_report_txt)
        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "", species_math)

        datahandling.save_sample_component(db_sample_component, sample_component_file)

        db_sample["properties"]["species_detection"] = db_sample_component["summary"]
        db_sample["properties"]["detected_species"] = db_sample_component["summary"]["name_classified_species_1"]
        db_sample = datahandling.datadump_template(db_sample, folder, "", set_sample_species)

        datahandling.save_sample(db_sample, sample_file)

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
