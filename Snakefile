#!/usr/bin/env python3

# TODO: should refactor a lot of this Snakefile into a more snakemake orientated solution utilizing wildcards
import re
import sys
import os
import datetime
import pandas
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "scripts"))
import datahandling

configfile: os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")

#Saving the config
component = "serumqc"

datahandling.save_yaml(config, "serumqc_config.yaml")

components = config["components"]
run_folder = config["run_folder"]
sample_sheet = config["sample_sheet"]
group = config["group"]
partition = config["partition"]
global_threads = config["threads"]
global_memory_in_GB = config["memory"]


onsuccess:
    print("Workflow complete")
    shell("touch serumqc_successfully_initialized_on_" + str(datetime.datetime.now()).replace(" ", "_"))


onerror:
    print("Workflow error")
    shell("touch serumqc_failed_to_initialized_on_" + str(datetime.datetime.now()).replace(" ", "_"))


rule all:
    input:
        component + "/" + component + "_complete"


rule setup:
    output:
        folder = directory(component)
    shell:
        "mkdir {output}"


rule_name = "export_conda_env"
rule export_conda_env:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
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
        component + "/conda_env.yaml"
    shell:
        "conda env export 1> {output} 2> {log.err_file}"


rule_name = "generate_git_hash"
rule generate_git_hash:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
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
        component + "/git_hash.txt"
    shell:
        "git --git-dir {workflow.basedir}/.git rev-parse HEAD 1> {output} 2> {log.err_file}"


rule initialize_components:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = component + "/log/" + rule_name + ".out.log",
        err_file = component + "/log/" + rule_name + ".err.log",
    benchmark:
        component + "/benchmarks/" + rule_name + ".benchmark"
    message:
        "Running step: {rule}"
    # Dynamic
    input:
        git_hash = rules.generate_git_hash.output,
        conda_env = rules.export_conda_env.output
    output:
        touch(component + "/initialize_components_complete"),
    run:
        git_hash = str(input.git_hash)
        conda_env = str(input.conda_env)

        sys.stdout.write("Started {}\n".format(rule_name))
        component_db = {}
        with open(git_hash, "r") as git_info:
            git_hash = git_info.readlines()[0].strip()
            component_db["git_hash"] = git_hash
        component_db["conda_env"] = datahandling.load_yaml(conda_env)
        component_db["config"] = config
        for component_name in components.split(","):
            component_db["name"] = component_name.strip()
            datahandling.save_component(component_db, component + "/" + component_name + ".yaml")
        sys.stdout.write("Done {}\n".format(rule_name))


rule_name = "initialize_samples_from_run_folder"
rule initialize_samples_from_run_folder:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
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
        run_folder = run_folder,
    output:
        touch(component + "/initialize_samples_from_run_folder")
    run:
        run_folder = str(input.run_folder)

        sys.stdout.write("Started {}\n".format(rule_name))
        for file in sorted(os.listdir(run_folder)):
            result = re.search(config["read_pattern"], file)
            if result and os.path.isfile(os.path.realpath(os.path.join(run_folder, file))):
                sample_name = result.group("sample_name")
                sample_config = sample_name + "/sample.yaml"
                sample_db = datahandling.load_sample(sample_config)
                sample_db["name"] = sample_name
                sample_db[result.group("paired_read_number")] = os.path.realpath(os.path.join(run_folder, file))
                # sample_db[result.group("paired_read_number") + "_md5sum"] = md5sum(os.path.realpath(os.path.join(run_folder, file)))
                datahandling.save_sample(sample_db, sample_config)
        sys.stdout.write("Done {}\n".format(rule_name))


rule_name = "check__provided_sample_info"
rule check__provided_sample_info:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
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
        rules.initialize_samples_from_run_folder.output,
        run_folder = run_folder,
    output:
        sample_sheet_tsv = component + "/sample_sheet.tsv",
    params:
        sample_sheet
    run:
        sample_sheet = str(params)
        corrected_sample_sheet_tsv = str(output.sample_sheet_tsv)

        sys.stdout.write("Started {}\n".format(rule_name))
        if not os.path.isfile(sample_sheet):
            df = pandas.DataFrame()
            df.to_csv(corrected_sample_sheet_tsv, sep='\t', index=False)
            return 0
        if sample_sheet.endswith(".xlsx"):
            df = pandas.read_excel(sample_sheet)
        else:  # assume it's a tsv
            df = pandas.read_table(sample_sheet)
        item_rename_dict = {}
        badly_named_samples = df[df["SampleID"].str.contains("^[a-zA-Z0-9\-_]+$") == False]  # samples which fail this have inappropriate characters
        for item in badly_named_samples["SampleID"].tolist():
            sys.stdout.write("Renaming '{}' to '{}'".format(item, re.sub(r'[^a-zA-Z0-9\-_]', "_", str(item).strip())))
        for item in df["SampleID"].tolist():
            item_rename_dict[item] = re.sub(r'[^a-zA-Z0-9\-_]', "_", str(item).strip())
        df["SampleID"] = df["SampleID"].map(item_rename_dict)
        duplicates = {}
        for item in df["SampleID"].tolist():
            if df["SampleID"].tolist().count(item) > 1:
                sys.stdout.write("Duplicate SampleID's exist {}".format(item))
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
        df.to_csv(corrected_sample_sheet_tsv, sep='\t', index=False)
        sys.stdout.write("Done {}\n".format(rule_name))


rule_name = "set_samples_from_sample_info"
rule set_samples_from_sample_info:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
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
        corrected_sample_sheet_tsv = rules.check__provided_sample_info.output,
    output:
        touch(component + "/set_samples_from_sample_info")
    params:
        sample_sheet
    run:
        corrected_sample_sheet_tsv = str(input.corrected_sample_sheet_tsv)

        sys.stdout.write("Started {}\n".format(rule_name))
        config = datahandling.load_config()
        # TODO: handle no sample sheet
        df = pandas.read_table(corrected_sample_sheet_tsv)
        for index, row in df.iterrows():
            sample_config = row["SampleID"] + "/sample.yaml"
            sample_db = datahandling.load_sample(sample_config)
            sample_db["sample_sheet"] = {}
            for column in df:
                column_name = column
                if column in config["samplesheet_column_mapping"]:
                    column_name = config["samplesheet_column_mapping"][column]
                sample_db["sample_sheet"][column_name] = row[column]
            datahandling.save_sample(sample_db, sample_config)
        sys.stdout.write("Done {}\n".format(rule_name))


rule_name = "add_components_to_samples"
rule add_components_to_samples:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
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
        rules.set_samples_from_sample_info.output,
        component = component,
        run_folder = run_folder,
    output:
        touch(component + "/add_components_to_samples"),
        touch(rules.all.input),
    params:
        sample_sheet
    run:
        run_folder = str(input.run_folder)
        serumqc_folder = str(input.component)

        sys.stdout.write("Started {}\n".format(rule_name))
        config = datahandling.load_config()
        for file in sorted(os.listdir(run_folder)):
            result = re.search(config["read_pattern"], file)
            if result and os.path.isfile(os.path.realpath(os.path.join(run_folder, file))):
                sample_name = result.group("sample_name")
                sample_config = sample_name + "/sample.yaml"
                sample_db = datahandling.load_sample(sample_config)
                sample_db["component_ids"] = sample_db.get("component_ids", [])
                for component in components:
                    # TODO: fix this
                    sample_db["component_ids"].append(datahandling.load_component(os.path.join(serumqc_folder, component + ".yaml")).get("_id",))
                sample_db[result.group("paired_read_number")] = os.path.realpath(os.path.join(run_folder, file))
                # sample_db[result.group("paired_read_number") + "_md5sum"] = md5sum(os.path.realpath(os.path.join(run_folder, file)))
                datahandling.save_sample(sample_db, sample_config)
        sys.stdout.write("Done {}\n".format(rule_name))

# rule_name = "initialize_run"
# rule initialize_run:
#     # Static
#     message:
#         "Running step:" + rule_name
#     threads:
#         global_threads
#     resources:
#         memory_in_GB = global_memory_in_GB
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
#         init_complete = component + "/initialize_components_complete",
#         run_folder = run_folder,
#     output:
#         samplesheet = "sample_sheet.tsv",
#         output = "run.yaml"
#     params:
#         samplesheet = sample_sheet,
#         partition = partition,
#         components = components,
#         group = group,
#         config = config,
#     script:
#         os.path.join(os.path.dirname(workflow.snakefile), "scripts/initialize_run.py")


# rule get_git_hash_of_serumqc:
#     input:
#         run_info_yaml_path = "run.yaml"
#     output:
#         git_hash = "serumqc/git_hash.txt"
#     run:
#         run_info = datahandling.load_run(input.run_info_yaml_path)
#         shell("git --git-dir {workflow.basedir}/.git rev-parse snakemake 1> {output}")
#         with open(output.git_hash, "r") as git_info:
#             git_hash = git_info.readlines()[0].strip()
#         run_info["run"]["git_hash"] = git_hash
#         datahandling.save_run(run_info, input.run_info_yaml_path)


# rule get_conda_env:
#     input:
#         git_hash = "serumqc/git_hash.txt",
#         run_info_yaml_path = "run.yaml"
#     output:
#         conda_yaml = "serumqc/conda.yaml"
#     run:
#         run_info = datahandling.load_run(input.run_info_yaml_path)
#         shell("conda env export 1> {output}")
#         run_info["run"]["conda_env"] = datahandling.load_yaml(output.conda_yaml)
#         datahandling.save_run(run_info, input.run_info_yaml_path)


# rule add_components_data_entry:
#     input:
#         git_hash = "serumqc/git_hash.txt",
#         conda_yaml = "serumqc/conda.yaml",
#     output:
#         components_db = ""

# rule create_end_file:
#     input:
#         "run.yaml"
#     output:
#         rules.all.input
#     shell:
#         """
#         bash run_cmd_serumqc.sh
#         touch {output}
#         """

# can break this down to 2 parts where you create the sample_sheet in one and then prep for run with the other
# rule start_run:
#     input:
#         "run.yaml"
#     output:
#         touch("run_started")
#     shell:
#         "bash run_cmd_serumqc.sh"