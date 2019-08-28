#!/usr/bin/env python3
"""
Initialization program for paired end Illumina reads
"""
# TODO: should refactor a lot of this Snakefile into a more snakemake orientated solution utilizing wildcards
import re
import sys
import os
import datetime
import pandas
import pkg_resources
import hashlib
import traceback
import shutil
from bifrostlib import datahandling
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "scripts"))

configfile: os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")

# hacky fix for AB#371, need to REDO
if config.get("samples_to_include", None) is not None:
    if type(config["samples_to_include"]) is tuple:
        config["samples_to_include"] = ",".join([str(i) for i in list(config["samples_to_include"])])
    if type(config["samples_to_include"]) is int:
        config["samples_to_include"] = str(config["samples_to_include"])

# Saving the config
component = "bifrost"
rerun_folder = component + "/delete_to_update"

datahandling.save_yaml(config, "config.yaml")
components = config["components"].split(",")
# raw_data_folder = config["raw_data_folder"]
# rename_samples = config["rename_samples"]
sample_folder = config["sample_folder"]
sample_sheet = config["sample_sheet"]
group = config["group"]
partition = config["partition"]
num_of_threads, memory_in_GB = config["threads"], config["memory"]


onsuccess:
    print("Workflow complete")
    shell("touch " + component + "_successfully_initialized_on_" + str(datetime.datetime.now()).replace(" ", "_"))


onerror:
    print("Workflow error")
    shell("touch " + component + "_failed_to_initialized_on_" + str(datetime.datetime.now()).replace(" ", "_"))


rule all:
    input:
        rerun_folder + "/" + component + "_complete"


rule setup:
    output:
        folder = directory(component)
    shell:
        "mkdir {output};"

rule setup_rerun:
    input:
        component
    output:
        rerun_folder = directory(component + "/delete_to_update")
    shell:
        "mkdir {output};"


rule make_components_dir:
    input:
        component
    output:
        folder = directory("components")
    shell:
        "mkdir {output};"

# TODO: temporarily shelved idea for anonymizing samples
# rule_name = "create_sample_folder"
# rule create_sample_folder:
#     # Static
#     message:
#         "Running step:" + rule_name
#     threads:
#         num_of_threads#     resources:
#         memory_in_GB
#     log:
#         out_file = component + "/log/" + rule_name + ".out.log",
#         err_file = component + "/log/" + rule_name + ".err.log",
#     benchmark:
#         component + "/benchmarks/" + rule_name + ".benchmark"
#     message:
#         "Running step: {rule}"
#     # Dynamic
#     input:
#         component,
#         raw_data_folder = raw_data_folder
#     output:
#         sample_folder = directory(sample_folder)
#     params:
#         rename_samples = rename_samples
#     run:
#         rename_samples = bool(params.rename_samples)
#         raw_data_folder = str(input.raw_data_folder)
#         sample_folder = str(output.sample_folder)
#         print(rename_samples, type(rename_samples))
#         if rename_samples is False:
#             shell("ln -s {raw_data_folder} {sample_folder}")
#         else:
#             shell("mkdir {sample_folder}")
#             i = 0
#             for file in sorted(os.listdir(raw_data_folder)):
#                 print(file)
#                 result = re.search(config["read_pattern"], file)
#                 if result and os.path.isfile(os.path.realpath(os.path.join(raw_data_folder, file))):
#                     i = i + 1
#                     new_sample_name = "SSI{}_R1.fastq.gz".format(i)
#                     print(new_sample_name)
#                     shell("ln -s {} {};".format( os.path.realpath(os.path.join(raw_data_folder, file)), os.path.join(sample_folder, new_sample_name)))


rule_name = "copy_run_info"
rule copy_run_info:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        component
    output:
        touch(rerun_folder + "/copy_run_info"),
        touch(component + "/copy_run_info_complete")
    params:
        sample_folder
    shell:
        """
        if [ -d \"{params}/InterOp\" ]; then cp -TR {params}/InterOp {input}/InterOp; fi;
        if [ -f \"{params}/RunInfo.xml\" ]; then cp {params}/RunInfo.xml {input}/RunInfo.xml; fi;
        if [ -f \"{params}/RunParams.xml\" ]; then cp {params}/RunParams.xml {input}/RunParams.xml; fi;
        """


rule_name = "initialize_components"
rule initialize_components:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        component = component,
        components_dir = rules.make_components_dir.output.folder
    output:
        touch(rerun_folder + "/initialize_components"),
        touch(component + "/initialize_components_complete"),
    params:
        rule_name = rule_name
    run:
        try:
            rule_name = str(params.rule_name)
            component = str(input.component)
            components_dir = str(input.components_dir)
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))

            for component_name in components:
                component_file = os.path.join(components_dir, component_name + ".yaml")
                if not os.path.isfile(component_file):
                    shutil.copyfile(os.path.join(os.path.dirname(workflow.snakefile), "components", component_name, "config.yaml"), component_file)
                db_component = datahandling.load_component(component_file)
                datahandling.save_component(db_component, component_file)

            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(traceback.format_exc()))


rule_name = "initialize_samples_from_sample_folder"
rule initialize_samples_from_sample_folder:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        component,
        sample_folder = sample_folder,
    output:
        touch(rerun_folder + "/initialize_samples_from_sample_folder"),
        touch(component + "/initialize_samples_from_sample_folder")
    params:
        rule_name = rule_name
    run:
        try:
            rule_name = str(params.rule_name)
            sample_folder = str(input.sample_folder)
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))

            unique_sample_names = {}
            for file in sorted(os.listdir(sample_folder)):
                result = re.search(config["read_pattern"], file)
                if result and os.path.isfile(os.path.realpath(os.path.join(sample_folder, file))):
                    sample_name = str(result.group("sample_name"))
                    unique_sample_names[sample_name] = unique_sample_names.get(sample_name, 0) + 1

            for sample_name in unique_sample_names:
                if config.get("samples_to_include", None) is None or sample_name in config["samples_to_include"].split(","):
                    sample_config = sample_name + "/sample.yaml"
                    sample_db = datahandling.load_sample(sample_config)
                    sample_db["name"] = sample_name
                    sample_db["reads"] = sample_db.get("reads", {})
                    sample_db["path"] = os.path.realpath(sample_name)
                    for file in sorted(os.listdir(sample_folder)):
                        result = re.search(config["read_pattern"], file)
                        if result and os.path.isfile(os.path.realpath(os.path.join(sample_folder, file))):
                            if sample_name == str(result.group("sample_name")):
                                sample_db["reads"]["R"+result.group("paired_read_number")] = os.path.realpath(os.path.join(sample_folder, file))
                                md5sum_key = result.group("paired_read_number") + "_md5sum"
                                if "md5skip" in config and config["md5skip"] and md5sum_key in sample_db["reads"]:
                                    pass
                                else:
                                    with open(os.path.realpath(os.path.join(sample_folder, file)), "rb") as fh:
                                        md5sum = hashlib.md5()
                                        for data in iter(lambda: fh.read(4096), b""):
                                            md5sum.update(data)
                                        sample_db["reads"][result.group("paired_read_number") + "_md5sum"] = md5sum.hexdigest()
                        sample_db["properties"] = sample_db.get("properties", {})
                        sample_db["report"] = sample_db.get("report", {})
                    datahandling.save_sample(sample_db, sample_config)
            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(traceback.format_exc()))


rule_name = "check__provided_sample_info"
rule check__provided_sample_info:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        rules.initialize_samples_from_sample_folder.output,
    output:
        touch(rerun_folder + "/check__provided_sample_info"),
        sample_sheet_tsv = component + "/sample_sheet.tsv",
    params:
        rule_name = rule_name,
        sample_sheet = sample_sheet
    run:
        try:
            rule_name = str(params.rule_name)
            sample_sheet = str(params.sample_sheet)
            corrected_sample_sheet_tsv = str(output.sample_sheet_tsv)
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            if not os.path.isfile(sample_sheet):
                df = pandas.DataFrame()
                df.to_csv(corrected_sample_sheet_tsv, sep="\t", index=False)
                return 0
            if sample_sheet.endswith(".xlsx"):
                df = pandas.read_excel(sample_sheet)
            else:  # assume it's a tsv
                df = pandas.read_table(sample_sheet)
            
            #Convert sample name to string in case all sample names are numbers
            df["SampleID"] = df["SampleID"].astype(str)
            # Dropping rows with no sample name.
            noname_index = df[df["SampleID"].isna()].index
            df.drop(noname_index, inplace=True)

            item_rename_dict = {}
            badly_named_samples = df[df["SampleID"].str.contains("^[a-zA-Z0-9\-_]+$") == False]  # samples which fail this have inappropriate characters
            for item in badly_named_samples["SampleID"].tolist():
                datahandling.log(log_out, "Renaming '{}' to '{}'".format(item, re.sub(r'[^a-zA-Z0-9\-_]', "_", str(item).strip())))
            for item in df["SampleID"].tolist():
                item_rename_dict[item] = re.sub(r'[^a-zA-Z0-9\-_]', "_", str(item).strip())
            df["SampleID"] = df["SampleID"].map(item_rename_dict)
            duplicates = {}
            for item in df["SampleID"].tolist():
                if df["SampleID"].tolist().count(item) > 1:
                    datahandling.log(log_out, "Duplicate SampleID's exist {}".format(item))
                    duplicates[item] = df["SampleID"].tolist().count(item)
                    # duplicates have to be done on id basis
            ridiculous_count = 0
            for i, row in df.iterrows():
                if df.loc[i, "SampleID"] in duplicates:
                    if df.loc[i, "SampleID"] + "_RENAMED-" + str(duplicates[df.loc[i, "SampleID"]]) in df["SampleID"].tolist():
                        ridiculous_count += 1
                        df.loc[i, "SampleID"] = "VERY_BAD_SAMPLE_NAME_" + str(ridiculous_count)
                    else:
                        df.loc[i, "SampleID"] = df.loc[i, "SampleID"] + "_RENAMED-" + str(duplicates[df.loc[i, "SampleID"]])
                    duplicates[row["SampleID"]] -= 1
            df.to_csv(corrected_sample_sheet_tsv, sep="\t", index=False)
            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(traceback.format_exc()))


rule_name = "set_samples_from_sample_info"
rule set_samples_from_sample_info:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        corrected_sample_sheet_tsv = rules.check__provided_sample_info.output.sample_sheet_tsv,
    output:
        touch(rerun_folder + "/set_samples_from_sample_info"),
        touch(component + "/set_samples_from_sample_info")
    params:
        rule_name = rule_name
    run:
        try:
            rule_name = str(params.rule_name)
            corrected_sample_sheet_tsv = str(input.corrected_sample_sheet_tsv)
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            config = datahandling.load_config()
            try:
                df = pandas.read_table(corrected_sample_sheet_tsv)
                #Convert sample name to string in case all sample names are numbers
                df["SampleID"] = df["SampleID"].astype(str)
                unnamed_sample_count = 0
                for index, row in df.iterrows():
                    sample_config = row["SampleID"] + "/sample.yaml"
                    # NOTE Sample name can be changed to Unnamed_ following the logic below (no name in sample sheet) however the .yaml path uses the value in the sample sheet
                    if config.get("samples_to_include", None) is None or row["SampleID"] in config["samples_to_include"].split(","):
                        sample_db = datahandling.load_sample(sample_config)
                        sample_db["sample_sheet"] = {}
                        for column in df:
                            column_name = column
                            for rename_column in config["samplesheet_column_mapping"]:
                                if config["samplesheet_column_mapping"][rename_column] == column:
                                    column_name = rename_column
                            sample_db["sample_sheet"][column_name] = row[column]
                        # If sample has no name (most likely because there were no reads
                        # in the sample folder) we have to specify a name.
                        if "name" not in sample_db:
                            if "sample_name" in sample_db["sample_sheet"]:
                                sample_db["name"] = sample_db["sample_sheet"]["sample_name"]
                            else:
                                # Increasing counter before so it starts with Unnamed_1
                                unnamed_sample_count += 1
                                sample_db["name"] = "Unnamed_" + unnamed_sample_count
                        sample_db["path"] = os.path.realpath(sample_db["name"])
                        datahandling.save_sample(sample_db, sample_config)
            except pandas.io.common.EmptyDataError:
                datahandling.log(log_err, ("No samplesheet data\n"))
            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(traceback.format_exc()))


rule_name = "set_sample_species"
rule set_sample_species:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        rules.set_samples_from_sample_info.output,
        corrected_sample_sheet_tsv = rules.check__provided_sample_info.output.sample_sheet_tsv,
    output:
        touch(rerun_folder + "/set_sample_species"),
        touch(component + "/set_sample_species")
    params:
        rule_name = rule_name
    run:
        try:
            rule_name = str(params.rule_name)
            corrected_sample_sheet_tsv = str(input.corrected_sample_sheet_tsv)
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            config = datahandling.load_config()
            try:
                df = pandas.read_table(corrected_sample_sheet_tsv)
                #Convert sample name to string in case all sample names are numbers
                df["SampleID"] = df["SampleID"].astype(str)
                for index, row in df.iterrows():
                    sample_config = row["SampleID"] + "/sample.yaml"
                    if config.get("samples_to_include", None) is None or row["SampleID"] in config["samples_to_include"].split(","):
                        sample_db = datahandling.load_sample(sample_config)

                        sample_db["properties"] = sample_db.get("properties", {})
                        provided_species = sample_db["sample_sheet"].get("provided_species")
                        if pandas.isna(provided_species):
                            provided_species = None
                        else:
                            species_db = datahandling.get_ncbi_species(
                                provided_species)
                            if species_db is None:
                                provided_species = str(provided_species)
                            else:
                                provided_species = species_db # Use proper name if exists.
                        sample_db["properties"]["provided_species"] = provided_species
                        datahandling.save_sample(sample_db, sample_config)

            except pandas.io.common.EmptyDataError:
                datahandling.log(log_err, "No samplesheet data\n")
            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(traceback.format_exc()))


rule_name = "add_components_to_samples"
rule add_components_to_samples:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        rules.initialize_components.output,
        rules.set_sample_species.output,
        component = component,
        sample_folder = sample_folder
    output:
        touch(rerun_folder + "/add_components_to_samples"),
        touch(component + "/add_components_to_samples"),
    params:
        rule_name = rule_name
    run:
        try:
            rule_name = str(params.rule_name)
            sample_folder = str(input.sample_folder)
            component = str(input.component)
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            config = datahandling.load_config()
            unique_sample_names = {}

            # I changed this to include all samples, even those with no data.
            # So that they can run components. (only stamps for now)
            for folder in sorted(os.listdir(".")):
                if os.path.isfile(os.path.realpath(os.path.join(folder, "sample.yaml"))):
                    sample_name = str(folder)
                    unique_sample_names[sample_name] = unique_sample_names.get(
                        sample_name, 0) + 1

            for sample_name in unique_sample_names:
                if config.get("samples_to_include", None) is None or sample_name in config["samples_to_include"].split(","):
                    sample_config = sample_name + "/sample.yaml"
                    sample_db = datahandling.load_sample(sample_config)
                    sample_db["components"] = sample_db.get("components", [])
                    for component_name in components:
                        component_id = datahandling.load_component("components/" + component_name + ".yaml").get("_id",)
                        if component_id is not None:
                            insert_component = True
                            for sample_component in sample_db["components"]:
                                if component_id == sample_component["_id"]:
                                    insert_component = False
                            if insert_component is True:
                                sample_db["components"].append({"name": component_name, "_id": component_id})
                    datahandling.save_sample(sample_db, sample_config)
            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(traceback.format_exc()))


rule_name = "initialize_sample_components_for_each_sample"
rule initialize_sample_components_for_each_sample:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        rules.add_components_to_samples.output,
        component = component,
        sample_folder = sample_folder
    output:
        touch(rerun_folder + "/initialize_sample_components_for_each_sample"),
        touch(component + "/initialize_sample_components_for_each_sample"),
    params:
        rule_name = rule_name
    run:
        try:
            rule_name = str(params.rule_name)
            sample_folder = str(input.sample_folder)
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            config = datahandling.load_config()
            unique_sample_names = {}
            # I changed this to include all samples, even those with no data.
            # So that they can run components. (only stamps for now)
            for folder in sorted(os.listdir(".")):
                if os.path.isfile(os.path.realpath(os.path.join(folder, "sample.yaml"))):
                    sample_name = str(folder)
                    unique_sample_names[sample_name] = unique_sample_names.get(
                        sample_name, 0) + 1

            for sample_name in unique_sample_names:
                if config.get("samples_to_include", None) is None or sample_name in config["samples_to_include"].split(","):
                    sample_config = sample_name + "/sample.yaml"
                    sample_db = datahandling.load_sample(sample_config)

                    sample_id = sample_db.get("_id",)
                    component_dict = sample_db.get("components",[])
                    for item in component_dict:
                        if item.get("name",) in components:
                            component_name = item.get("name",)
                            component_id = item.get("_id",)
                            sample_component_path = sample_name + "/" + sample_name + "__" + component_name + ".yaml"
                            sample_component_folder_path = os.path.realpath(
                                os.path.join(sample_name, component_name))
                            sample_component_db = datahandling.load_sample_component(sample_component_path)
                            sample_component_db["sample"] = {"name": sample_name, "_id": sample_id}
                            sample_component_db["component"] = {"name": component_name, "_id": component_id}
                            sample_component_db["status"] = "initialized"
                            sample_component_db["path"] = sample_component_folder_path
                            datahandling.save_sample_component(sample_component_db, sample_component_path)
            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(traceback.format_exc()))


rule_name = "initialize_run"
rule initialize_run:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        rules.copy_run_info.output,
        rules.initialize_sample_components_for_each_sample.output,
        sample_folder = sample_folder,
        run_folder = component
    output:
        touch(rerun_folder + "/initialize_run"),
        touch(component + "/initialize_run"),
    params:
        rule_name = rule_name
    run:
        try:
            rule_name = str(params.rule_name)
            sample_folder = str(input.sample_folder)
            run_folder = str(input.run_folder)
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            config = datahandling.load_config()
            unique_sample_names = {}

            run_db = datahandling.load_run(run_folder + "/run.yaml")
            run_db["name"] = config.get("run_name", os.getcwd().split("/")[-1])
            if "type" in config:
                run_db["type"] = config["type"]
            else:
                if "type" not in run_db:
                    run_db["type"] = "default"
            run_db["path"] = os.path.realpath(run_folder)
            for folder in sorted(os.listdir(".")):
                if os.path.isfile(os.path.realpath(os.path.join(folder, "sample.yaml"))):
                    sample_name = str(folder)
                    unique_sample_names[sample_name] = unique_sample_names.get(sample_name, 0) + 1

            run_db["samples"] = run_db.get("samples", [])
            # Todo: handle change in samples properly
            for sample_name in unique_sample_names:
                if config.get("samples_to_include", None) is None or sample_name in config["samples_to_include"].split(","):
                    sample_config = sample_name + "/sample.yaml"
                    sample_db = datahandling.load_sample(sample_config)
                    sample_id = sample_db.get("_id",)
                    if sample_id is not None:
                        insert_sample = True
                        for sample in run_db["samples"]:
                            if sample_id == sample["_id"]:
                                insert_sample = False
                        if insert_sample is True:
                            run_db["samples"].append({"name": sample_name, "_id": sample_id})

            run_db["components"] = run_db.get("components", [])
            for component_name in components:
                component_id = datahandling.load_component("components/" + component_name + ".yaml").get("_id",)
                if component_id is not None:
                    insert_component = True
                    for run_components in run_db["components"]:
                        if component_id == run_components["_id"]:
                            insert_component = False
                    if insert_component is True:
                        run_db["components"].append({"name": component_name, "_id": component_id})

            datahandling.save_run(run_db, component + "/run.yaml")
            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(traceback.format_exc()))


rule_name = "setup_sample_components_to_run"
rule setup_sample_components_to_run:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        num_of_threads
    resources:
        memory_in_GB = memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        rules.initialize_run.output,
        component = component,
        sample_folder = sample_folder
    output:
        touch(rerun_folder + "/setup_sample_components_to_run"),
        bash_file = "run_cmd_" + component + ".sh"
    params:
        rule_name = rule_name
    run:
        try:
            rule_name = str(params.rule_name)
            sample_folder = str(input.sample_folder)
            component = str(input.component)
            run_cmd = str(output.bash_file)
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            config = datahandling.load_config()
            unique_sample_names = {}
            # I changed this to include all samples, even those with no data.
            # So that they can run components. (only stamps for now)
            for folder in sorted(os.listdir(".")):
                if os.path.isfile(os.path.realpath(os.path.join(folder, "sample.yaml"))):
                    sample_name = str(folder)
                    unique_sample_names[sample_name] = unique_sample_names.get(
                        sample_name, 0) + 1

            with open(run_cmd, "w") as run_cmd_handle:
                run_cmd_handle.write("bifrost__job_ids=;")
                for sample_name in unique_sample_names:
                    if config.get("samples_to_include", None) is None or sample_name in config["samples_to_include"].split(","):
                        current_time = datetime.datetime.now()
                        with open(sample_name + "/cmd_" + component + "_{}.sh".format(current_time), "w") as command:
                            sample_config = sample_name + "/sample.yaml"
                            sample_db = datahandling.load_sample(sample_config)

                            # Default priority
                            partition = config["partition"]

                            # Set sample priority to value from run_metadata then map it based off config file
                            if "sample_sheet" in sample_db:
                                if "priority" in sample_db["sample_sheet"]:
                                    if isinstance(sample_db["sample_sheet"]["priority"], str) is True and sample_db["sample_sheet"]["priority"].lower() in config["priority_mapping"]:
                                        partition = config["priority_mapping"][sample_db["sample_sheet"]["priority"].lower()]

                            # Overrule run_metadata and default value
                            if "force_partition" in config:
                                partition = config["force_partition"]

                            command.write("#!/bin/sh\n")
                            if config["grid"] == "torque":
                                if "advres" in config and config["advres"]:
                                    advres = ",advres={}".format(config["advres"])
                                else:
                                    advres = ''
                                torque_node = ",nodes=1:ppn={}".format(config["threads"])

                                command.write("#PBS -V -d . -w . -l mem={}gb{},walltime={}{} -N '{}_{}' -W group_list={} -A {} \n".format(config["memory"], torque_node, config["walltime"], advres, component, sample_name, group, group))
                            elif config["grid"] == "slurm":
                                command.write("#SBATCH --mem={}G -p {} -c {} -t {} -J '{}_{}'\n".format(
                                    config["memory"], partition, config["threads"], config["walltime"], component, sample_name))

                            if "tmp_dir" in config:
                                tmp_dir = " --shadow-prefix {}".format(
                                    config["tmp_dir"])
                            else:
                                tmp_dir = ""

                            for component_name in components:
                                component_file = os.path.dirname(workflow.snakefile) + "/components/" + component_name + "/pipeline.smk"
                                if sample_name in config["samples_to_ignore"]:
                                    sample_component_db = datahandling.load_sample_component(sample_name + "/" + sample_name + "__" + component_name + ".yaml")
                                    sample_component_db["status"] = "skipped"
                                    sample_component_db["setup_date"] = current_time
                                    datahandling.save_sample_component(sample_component_db, sample_name + "/" + sample_name + "__" + component_name + ".yaml")
                                elif os.path.isfile(component_file):
                                    unlock = ""
                                    if config["unlock"]:
                                        unlock = "--unlock"
                                    else:
                                        # Only delete directory on non unlock mode
                                        command.write(
                                            "if [ -d \"{}\" ]; then rm -r {}; fi;\n".format(component_name, component_name))
                                    sample_component_db = datahandling.load_sample_component(sample_name + "/" + sample_name + "__" + component_name + ".yaml")
                                    sample_component_db["status"] = "queued to run"
                                    sample_component_db["setup_date"] = current_time
                                    datahandling.save_sample_component(sample_component_db, sample_name + "/" + sample_name + "__" + component_name + ".yaml")
                                    command.write("snakemake {} --use-singularity  --singularity-args \"{}\" --singularity-prefix \"{}\" --restart-times {} --cores {} -s {} {} --config Sample={} sample_id={} component_id={}; \n".format(tmp_dir, "-B " + sample_db["reads"]["R1"] + "," + sample_db["reads"]["R2"] + "," + os.getcwd() + "," + os.path.dirname(os.getenv("BIFROST_DB_KEY")), config["singularity_prefix"], config["restart_times"], config["threads"], component_file, unlock, "sample.yaml", sample_component_db["sample"]["_id"], sample_component_db["component"]["_id"]))
                                else:
                                    datahandling.log(log_err, "Error component not found:{} {}".format(component_name, component_file))
                                    sample_component_db = datahandling.load_sample_component(sample_name + "/" + sample_name + "__" + component_name + ".yaml")
                                    sample_component_db["status"] = "component missing"
                                    sample_component_db["setup_date"] = current_time
                                    datahandling.save_sample_component(sample_component_db, sample_name + "/" + sample_name + "__" + component_name + ".yaml")

                        os.chmod(os.path.join(sample_name, "cmd_{}_{}.sh".format(component, current_time)), 0o777)
                        if os.path.islink(os.path.join(sample_name, "cmd_{}.sh").format(component)):
                            os.remove(os.path.join(sample_name, "cmd_{}.sh").format(component))
                        os.symlink(os.path.realpath(os.path.join(sample_name, "cmd_{}_{}.sh".format(component, current_time))), os.path.join(sample_name, "cmd_" + component + ".sh"))
                        run_cmd_handle.write("cd {};\n".format(sample_name))
                        if config["grid"] == "torque":
                            run_cmd_handle.write("[ -z \"$bifrost__job_ids\" ] && bifrost__job_ids=$(qsub -h cmd_{}.sh) || bifrost__job_ids=$bifrost__job_ids:$(qsub -h cmd_{}.sh);\n".format(component, component))  # dependent on grid engine
                        elif config["grid"] == "slurm":
                            run_cmd_handle.write("[ -z \"$bifrost__job_ids\" ] && bifrost__job_ids=$(sbatch --hold --parsable cmd_{}.sh) || bifrost__job_ids=$bifrost__job_ids:$(sbatch --hold --parsable cmd_{}.sh);\n".format(component, component))  # dependent on grid engine
                        else:
                            run_cmd_handle.write("bash cmd_{}.sh;\n".format(component))
                        run_cmd_handle.write("cd {};\n".format(os.getcwd()))
                if os.path.isfile(os.path.join(os.path.dirname(workflow.snakefile), "scripts/final_script.sh")):
                    if config["grid"] == "torque":
                        run_cmd_handle.write("[ ! -z \"$bifrost__job_ids\" ] && bifrost__job_ids=$bifrost__job_ids:$(qsub -h -w depend=afterany:$bifrost__job_ids {}.sh);\n".format(os.path.join(os.path.dirname(workflow.snakefile), "scripts/final_script.sh")))  # dependent on grid engine
                    elif config["grid"] == "slurm":
                        run_cmd_handle.write("[ ! -z \"$bifrost__job_ids\" ] && bifrost__job_ids=$bifrost__job_ids:$(sbatch --hold -p {} --parsable -d afterany:$bifrost__job_ids {});\n".format(partition, os.path.join(os.path.dirname(workflow.snakefile), "scripts/final_script.sh")))  # dependent on grid engine
                    else:
                        run_cmd_handle.write("bash {};\n".format(os.path.join(os.path.dirname(workflow.snakefile), "scripts/final_script.sh")))
                if config["grid"] == "torque":
                    run_cmd_handle.write("qrls ${bifrost__job_ids//:/ };\n")
                elif config["grid"] == "slurm":
                    run_cmd_handle.write("scontrol release JobId=${bifrost__job_ids//:/ };\n")

            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(traceback.format_exc()))


rule create_end_file:
    input:
        rules.setup_sample_components_to_run.output,
        rerun_folder + "/initialize_samples_from_sample_folder"
    output:
        rules.all.input
    params:
        bashcmd = "run_cmd_{}.sh".format(component)
    run:
        bashcmd = params.bashcmd
        output = output
        if not config["init_only"]:
            shell("bash {} && touch {}".format(bashcmd, output))
        else:
            shell("touch {}".format(output))
