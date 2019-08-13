import pkg_resources
import datetime
import os
import re
import sys
import traceback
from bifrostlib import datahandling

config = datahandling.load_config()


def extract_cge_resfinder_data(file_path, key, data_dict):
    buffer = datahandling.load_yaml(file_path)
    data_dict["results"] = buffer
    return data_dict


def convert_summary_for_reporter(data_dict):
    for anti_biotic_class in data_dict["results"]["resfinder"]["results"]:
        for gene in anti_biotic_class:
            data_dict["reporter"]["content"].append([gene["resistance_gene"], gene["coverage"], gene["identity"], anti_biotic_class, gene["predicted_phenotype"]])  # table rows
    return data_dict


def script__datadump(output, folder, sample_file, component_file, sample_component_file, log):
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
        db_sample_component["reporter"] = db_component["db_values_changes"]["sample"]["reporter"]["resistance"]

        # Data extractions
        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "data_resfinder.json", extract_cge_resfinder_data)
        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "", convert_summary_for_reporter)

        # Save to sample component
        datahandling.save_sample_component(db_sample_component, sample_component_file)
        # Save summary and reporter results into sample
        db_sample["properties"]["resistance"] = db_sample_component["summary"]
        db_sample["reporter"]["resistance"] = db_sample_component["reporter"]
        datahandling.save_sample(db_sample, sample_file)

        open(output, 'w').close()  # touch file

    except Exception:
        datahandling.log(log_out, "Exception in {}\n".format(this_function_name))
        datahandling.log(log_err, str(traceback.format_exc()))
        raise Exception
        return 1

    finally:
        datahandling.log(log_out, "Done {}\n".format(this_function_name))
        return 0


script__datadump(
    snakemake.output.complete,
    snakemake.params.folder,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.sample_component_file,
    snakemake.log)
