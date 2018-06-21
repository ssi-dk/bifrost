#!/usr/bin/env python3

# TODO: should refactor a lot of this Snakefile into a more snakemake orientated solution utilizing wildcards

import re
import sys
import os
import datetime

configfile: os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")

components = str(config["components"])
run_folder = str(config["run_folder"])
sample_sheet = str(config["sample_sheet"])
group = str(config["group"])

partition = str(config["partition"])
global_threads = config["threads"]
global_memory_in_GB = config["memory"]
folder_name = "run_info"
# my understanding is all helps specify final output
onsuccess:
    print("Workflow complete")
    shell("touch qcquickie_successfully_initialized_on_" + str(datetime.datetime.now()).replace(" ", "_"))

onerror:
    print("Workflow error")
    shell("touch qcquickie_failed_to_initialized_on_" + str(datetime.datetime.now()).replace(" ", "_"))

ruleorder: setup > initialize_run


rule all:
    input:
        "run.yaml"


rule setup:
    output:
        directory = folder_name
    shell:
        "mkdir {output}"


rule initialize_run:
    message:
        "Running step: {rule}"
    input:
        run_folder = run_folder,
    output:
        samplesheet = "sample_sheet.tsv",
        output = "run.yaml"
    params:
        samplesheet = sample_sheet,
        partition = partition,
        components = components,
        group = group,
        config = config,
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    conda:
        "envs/python_packages.yaml"
    log:
        os.path.join(folder_name, "log/initialize_run.log")
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/initialize_run.py")


# rule setup_run_commands:
#     input:
#         sample_sheet = "sample_sheet.tsv"
#     output:
#         "run.yaml"

# rule check_provided_sample_sheet:
#     input:
#         directory
#     output:
#         sample_sheet = "sample_sheet.tsv"
#     params:
#         sample_sheet = "sample_sheet.xlsx"
#     run:
#         sys.stdout.write("Started check__provided_sample_info\n")
#         # change logic to apply to N_WGS and then used N_WGS to replace H_WGS
#         if not os.path.isfile(params.sample_sheet):
#             df = pandas.DataFrame()
#             df.to_csv(output.sample_sheet, sep='\t', index=False)
#             return 0
#         if sample_sheet.endswith(".xlsx"):
#             df = pandas.read_excel(params.sample_sheet)
#         else:  # assume it's a tsv
#             df = pandas.read_table(params.sample_sheet)
#         item_rename_dict = {}
#         badly_named_samples = df[df["SampleID"].str.contains("^[a-zA-Z0-9\-_]+$") == False]  # samples which fail this have inappropriate characters
#         for item in badly_named_samples["SampleID"].tolist():
#             print("Renaming '{}' to '{}'".format(item, re.sub(r'[^a-zA-Z0-9\-_]', "_", str(item).strip())))
#         for item in df["SampleID"].tolist():
#             item_rename_dict[item] = re.sub(r'[^a-zA-Z0-9\-_]', "_", str(item).strip())
#         df["SampleID"] = df["SampleID"].map(item_rename_dict)
#         duplicates = {}
#         for item in df["SampleID"].tolist():
#             if df["SampleID"].tolist().count(item) > 1:
#                 print("Duplicate SampleID's exist {}".format(item))
#                 duplicates[item] = df["SampleID"].tolist().count(item)
#                 # duplicates have to be done on id basis
#         ridiculous_count = 0
#         for i, row in df.iterrows():
#             if df.loc[i, "SampleID"] in duplicates:
#                 if df.loc[i, "SampleID"] + "_RENAMED-" + str(duplicates[df.loc[i, "SampleID"]]) in df["SampleID"].tolist():
#                     ridiculous_count += 1
#                     df.loc[i, "SampleID"] = "VERY_BAD_SAMPLE_NAME_" + str(ridiculous_count)
#                 else:
#                     df.loc[i, "SampleID"] = df.loc[i, "SampleID"] + "_RENAMED-" + str(duplicates[df.loc[i, "SampleID"]])
#                 duplicates[row["SampleID"]] -= 1
#         df.to_csv(output.sample_sheet, sep='\t', index=False)
#         sys.stdout.write("Done check__provided_sample_info\n")

# rule 
# rule 