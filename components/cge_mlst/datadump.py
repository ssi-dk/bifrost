import pkg_resources
import datetime
import os
import re
from bifrostlib import datahandling
import sys

config = datahandling.load_config()


def script__datadump(folder, sample_file, component_file, sample_component_file, log):
    db_sample = datahandling.load_sample(sample_file)
    db_component = datahandling.load_component(component_file)
    db_sample_component = datahandling.load_sample_component(sample_component_file)

    species = db_sample["properties"]["species"]

    db_sample_component["summary"] = {"db": [], "strain": [], "alleles": [], "component": {"id": db_component["_id"], "date": datetime.datetime.utcnow()}}

    mlst_species = db_component["mlst_species_mapping"][species]
    for mlst_entry in mlst_species:
        mlst_entry_db = datahandling.load_yaml(folder + "/" + mlst_entry + "/data.json")
        db_sample_component["results"][mlst_entry] = mlst_entry_db
        db_sample_component["summary"]["db"].append(mlst_entry)
        db_sample_component["summary"]["strain"].append(mlst_entry_db["mlst"]["results"].get("sequence_type", "NA"))
        db_sample_component["summary"]["alleles"].append(",".join([mlst_entry_db["mlst"]["results"]["allele_profile"][i]["allele_name"] for i in [i for i in mlst_entry_db["mlst"]["results"]["allele_profile"]]]))

    db_sample["properties"]["mlst"] = db_sample_component["summary"]
    datahandling.save_sample_component(db_sample_component, sample_component_file)
    datahandling.save_sample(db_sample, sample_file)

    return 0


script__datadump(
    snakemake.params.folder,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.sample_component_file,
    snakemake.log)
