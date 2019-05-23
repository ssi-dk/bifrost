import pkg_resources
import os
import re
from bifrostlib import datahandling
import sys

config = datahandling.load_config()

def script__datadump(folder, sample, sample_file_name, component_file_name):
    db_sample = datahandling.load_sample(sample_file_name)
    db_component = datahandling.load_component(component_file_name)

    folder = str(folder)
    sample = str(sample)

    datadump_dict = datahandling.load_sample_component(sample)
    datadump_dict["summary"] = datadump_dict.get("summary", {})
    datadump_dict["results"] = datadump_dict.get("results", {})

    species = db_sample["properties"]["species"]

    datadump_dict["summary"]["db"] = []
    datadump_dict["summary"]["strain"] = []
    datadump_dict["summary"]["alleles"] = []
    datadump_dict["summary"]["component"] = {"id" :db_component["_id"]}

    mlst_species = db_component["mlst_species_mapping"][species]
    for mlst_entry in mlst_species:
        mlst_entry_db = datahandling.load_yaml("cge_mlst/" + mlst_entry + "/data.json")
        datadump_dict["results"][mlst_entry] = mlst_entry_db
        datadump_dict["summary"]["db"].append(mlst_entry)
        datadump_dict["summary"]["strain"].append(mlst_entry_db["mlst"]["results"].get("sequence_type","NA"))
        datadump_dict["summary"]["alleles"].append(",".join([mlst_entry_db["mlst"]["results"]["allele_profile"][i]["allele_name"] for i in [i for i in mlst_entry_db["mlst"]["results"]["allele_profile"]]]))

    db_sample["properties"]["mlst"] = datadump_dict["summary"]
    datahandling.save_sample_component(datadump_dict, sample)
    datahandling.save_sample(db_sample, sample_file_name)

    return 0


script__datadump(snakemake.params.folder, snakemake.params.sample, snakemake.params.sample_file_name, snakemake.params.component_file_name)