import pkg_resources
import os
import re
from bifrostlib import datahandling
import sys

config = datahandling.load_config()


def extract_tsv(datadump_dict, folder, relative_path):
    relative_path_key = relative_path.replace(".", "_")
    if os.path.isfile(os.path.join(folder, relative_path)):
        datadump_dict["results"][relative_path_key] = {}
        try:
            with open(os.path.join(folder, relative_path), "r") as tsv_file:
                values = []
                headers = None
                for line in tsv_file:
                    line = line.strip()
                    if headers is None:
                        headers = line.split('\t')
                    else:
                        row = line.split('\t')
                        values.append(dict(zip(headers, row)))
        except Exception as e:
            sys.stderr.write(relative_path, e)
            sys.stderr.write(datadump_dict)
            datadump_dict["results"][relative_path_key]["status"] = "datadumper error"
        datadump_dict["results"][relative_path_key]["values"] = values
    return datadump_dict
    

def script__datadump_virulencefinder(folder, sample):
    folder = str(folder)
    sample = str(sample)

    datadump_dict = datahandling.load_sample_component(sample)
    datadump_dict["summary"] = datadump_dict.get("summary", {})
    datadump_dict["results"] = datadump_dict.get("results", {})

    datadump_dict = extract_tsv(
        datadump_dict, folder, "ariba_virulencefinder/report.tsv")

    datadump_dict = extract_tsv(
        datadump_dict, folder, "abricate_on_virulencefinder_from_ariba.tsv")
    
    # Summary:
    try:
        datadump_dict["summary"]["ariba_virulencefinder"] = datadump_dict["results"]["abricate_on_virulencefinder_from_ariba_tsv"]["values"]
    except KeyError as e:
        datadump_dict["summary"]["ariba_virulencefinder"] = "KeyError: {}".format(e)

    datahandling.save_sample_component(datadump_dict, sample)
    
    return 0

script__datadump_virulencefinder(snakemake.params.folder, snakemake.params.sample)
