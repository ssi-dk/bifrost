import pkg_resources
import datetime
import os
import re
import sys
from bifrostlib import datahandling

config = datahandling.load_config()


def extract_cge_resfinder_data(file_path, key, data_dict):
    buffer = datahandling.load_yaml(file_path)
    data_dict["results"] = buffer
    data_dict["summary"]["genes"].append(buffer["resfinder"]["results"].get("sequence_type", "NA"))
    return data_dict


def script__datadump(folder, sample_file, component_file, sample_component_file, log):
    db_sample = datahandling.load_sample(sample_file)
    db_component = datahandling.load_component(component_file)
    db_sample_component = datahandling.load_sample_component(sample_component_file)

    # Initialization of values to save into sample summary
    db_sample_component["summary"] = {"genes": [], "component": {"id": db_component["_id"], "date": datetime.datetime.utcnow()}}

    # Data extractions
    db_sample_component = datahandling.datadump_template(db_sample_component, folder, "/data.json", extract_cge_resfinder_data)

    # Save summarized results into sample
    db_sample["properties"]["resfinder"] = db_sample_component["summary"]
    datahandling.save_sample(db_sample, sample_file)
    datahandling.save_sample_component(db_sample_component, sample_component_file)

    return 0


script__datadump(
    snakemake.params.folder,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.sample_component_file,
    snakemake.log)
