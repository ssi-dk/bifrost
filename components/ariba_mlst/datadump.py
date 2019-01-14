import pkg_resources
import os
import re
from bifrostlib import datahandling
import sys

config = datahandling.load_config()


def extract_tsv(datadump_dict, analyzer_folder, relative_path):
    relative_path_key = relative_path.replace(".", "_")
    if os.path.isfile(os.path.join(analyzer_folder, relative_path)):
        datadump_dict["results"][relative_path_key] = {}
        try:
            with open(os.path.join(analyzer_folder, relative_path), "r") as tsv_file:
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


def script__datadump_ariba_mlst(analyzer_folder, sample):
    analyzer_folder = str(analyzer_folder)
    sample = str(sample)

    datadump_dict = datahandling.load_sample_component(sample)
    datadump_dict["summary"] = datadump_dict.get("summary", {})
    datadump_dict["results"] = datadump_dict.get("results", {})

    datadump_dict = extract_tsv(
        datadump_dict, analyzer_folder, "ariba_mlst/report.tsv")

    datadump_dict = extract_tsv(
        datadump_dict, analyzer_folder, "ariba_mlst/mlst_report.tsv")

    # Summary:
    try:
        datadump_dict["summary"]["mlst_report"] = ",".join(
            ["{}:{}".format(key, val) for key, val in
                datadump_dict["results"]["ariba_mlst/mlst_report_tsv"]["values"][0].items()])
    except KeyError as e:
        datadump_dict["summary"]["mlst_report"] = "KeyError: {}".format(e)

    datahandling.save_sample_component(datadump_dict, sample)

    return 0

script__datadump_ariba_mlst(snakemake.params.folder, snakemake.params.sample)
