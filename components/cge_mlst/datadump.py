import pkg_resources
import os
import re
from bifrostlib import datahandling
import sys

config = datahandling.load_config()

def script__datadump(folder, sample, sample_yaml):
    folder = str(folder)
    sample = str(sample)

    # datadump_dict = datahandling.load_sample_component(sample)
    # datadump_dict["summary"] = datadump_dict.get("summary", {})
    # datadump_dict["results"] = datadump_dict.get("results", {})
    # mlst_database = datahandling.get_mlst_species_DB(sample_yaml)
    # datadump_dict["results"]["mlst_db"] = mlst_database
    # datadump_dict["summary"]["mlst_db"] = mlst_database

    # datadump_dict = extract_tsv(
    #     datadump_dict, folder, "ariba_mlst/report.tsv")

    # datadump_dict = extract_tsv(
    #     datadump_dict, folder, "ariba_mlst/mlst_report.tsv")

    # # Summary:
    # try:
    #     datadump_dict["summary"]["mlst_report"] = ",".join(
    #         ["{}:{}".format(key, val) for key, val in
    #             datadump_dict["results"]["ariba_mlst/mlst_report_tsv"]["values"][0].items()])
    # except KeyError as e:
    #     datadump_dict["summary"]["mlst_report"] = "KeyError: {}".format(e)

    # datahandling.save_sample_component(datadump_dict, sample)

    return 0

script__datadump(snakemake.params.folder, snakemake.params.sample, snakemake.params.sample_yaml)
