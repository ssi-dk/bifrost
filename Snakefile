#!/usr/bin/env python3

# TODO: should refactor a lot of this Snakefile into a more snakemake orientated solution utilizing wildcards
import re
import sys
import os
import datetime
import pandas
import pkg_resources
import hashlib
import glob
sys.path.append(os.path.join(os.path.dirname(workflow.snakefile), "scripts"))
import datahandling

configfile: os.path.join(os.path.dirname(workflow.snakefile), "config.yaml")
#Saving the config
component = "serumqc"
rerun_folder = component + "/delete_to_update"

datahandling.save_yaml(config, "serumqc_config.yaml")

components = config["components"].split(",")
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
        rerun_folder + "/" + component + "_complete"


rule setup:
    output:
        folder = directory(component)
    shell:
        "mkdir {output}"

rule setup_rerun:
    input:
        component
    output:
        rerun_folder = directory(component + "/delete_to_update")
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
        touch(rerun_folder + "/export_conda_env"),
        conda_yaml = component + "/conda_env.yaml",
    shell:
        "conda env export 1> {output.conda_yaml} 2> {log.err_file}"


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
        touch(rerun_folder + "/generate_git_hash"),
        git_hash = component + "/git_hash.txt"
    shell:
        "git --git-dir {workflow.basedir}/.git rev-parse HEAD 1> {output.git_hash} 2> {log.err_file}"


rule_name = "copy_run_info"
rule copy_run_info:
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
        touch(rerun_folder + "/copy_run_info"),
        touch(component + "/copy_run_info_complete")
    params:
        run_folder
    shell:
        """
        if [ -d \"{params}/InterOp\" ]; then cp -TR {params}/InterOp {input}/InterOp ; fi;
        if [ -f \"{params}/RunInfo.xml\" ]; then cp {params}/RunInfo.xml {input}/RunInfo.xml; fi;
        if [ -f \"{params}/RunParams.xml\" ]; then cp {params}/RunParams.xml {input}/RunParams.xml; fi;
        """


rule_name = "initialize_components"
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
        component = component,
        git_hash = rules.generate_git_hash.output.git_hash,
        conda_env = rules.export_conda_env.output.conda_yaml,
    output:
        touch(rerun_folder + "/initialize_components"),
        touch(component + "/initialize_components_complete"),
    params:
        rule_name = rule_name
    run:
        rule_name = str(params.rule_name)
        git_hash = str(input.git_hash)
        conda_env = str(input.conda_env)
        component = str(input.component)

        sys.stdout.write("Started {}\n".format(rule_name))
        component_db = {}
        with open(git_hash, "r") as git_info:
            git_hash = git_info.readlines()[0].strip()
            component_db["git_hash"] = git_hash
        component_db["conda_env"] = datahandling.load_yaml(conda_env)
        component_db["config"] = config
        for component_name in components:
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
        touch(rerun_folder + "/initialize_samples_from_run_folder"),
        touch(component + "/initialize_samples_from_run_folder")
    params:
        rule_name = rule_name
    run:
        rule_name = str(params.rule_name)
        run_folder = str(input.run_folder)

        sys.stdout.write("Started {}\n".format(rule_name))

        unique_sample_names = {}
        for file in sorted(os.listdir(run_folder)):
            result = re.search(config["read_pattern"], file)
            if result and os.path.isfile(os.path.realpath(os.path.join(run_folder, file))):
                sample_name = result.group("sample_name")
                unique_sample_names[sample_name] = unique_sample_names.get(sample_name, 0) + 1

        for sample_name in unique_sample_names:
            sample_config = sample_name + "/sample.yaml"
            sample_db = datahandling.load_sample(sample_config)
            sample_db["name"] = sample_name
            sample_db["reads"] = sample_db.get("reads", {})
            files = glob.glob(os.path.realpath(os.path.join(run_folder, sample_name)) + "*")
            for file in files:
                result = re.search(config["read_pattern"], file)
                if result and os.path.isfile(os.path.realpath(os.path.join(run_folder, file))):
                    sample_db["reads"][result.group("paired_read_number")] = os.path.realpath(os.path.join(run_folder, file))
                # may be better to move this out
                    with open(os.path.realpath(os.path.join(run_folder, file)), "rb") as fh:
                        md5sum = hashlib.md5()
                        for data in iter(lambda: fh.read(4096), b""):
                            md5sum.update(data)
                        sample_db["reads"][result.group("paired_read_number") + "_md5sum"] = md5sum.hexdigest()
                sample_db["properties"] = {}  # init for others
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
        rules.initialize_samples_from_run_folder.output,
    output:
        touch(rerun_folder + "/check__provided_sample_info"),
        sample_sheet_tsv = component + "/sample_sheet.tsv",
    params:
        rule_name = rule_name,
        sample_sheet = sample_sheet
    run:
        rule_name = str(params.rule_name)
        sample_sheet = str(params.sample_sheet)
        corrected_sample_sheet_tsv = str(output.sample_sheet_tsv)

        sys.stdout.write("Started {}\n".format(rule_name))
        if not os.path.isfile(sample_sheet):
            df = pandas.DataFrame()
            df.to_csv(corrected_sample_sheet_tsv, sep="\t", index=False)
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
        df.to_csv(corrected_sample_sheet_tsv, sep="\t", index=False)
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
        corrected_sample_sheet_tsv = rules.check__provided_sample_info.output.sample_sheet_tsv,
    output:
        touch(rerun_folder + "/set_samples_from_sample_info"),
        touch(component + "/set_samples_from_sample_info")
    params:
        rule_name = rule_name
    run:
        rule_name = str(params.rule_name)
        corrected_sample_sheet_tsv = str(input.corrected_sample_sheet_tsv)

        sys.stdout.write("Started {}\n".format(rule_name))
        config = datahandling.load_config()
        # TODO: handle no sample sheet
        try:
            df = pandas.read_table(corrected_sample_sheet_tsv)
            unnamed_sample_count = 0
            for index, row in df.iterrows():
                sample_config = row["SampleID"] + "/sample.yaml"
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
                datahandling.save_sample(sample_db, sample_config)
        except pandas.io.common.EmptyDataError:
            sys.stderr.write("No samplesheet data\n")
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
        run_folder = run_folder
    output:
        touch(rerun_folder + "/add_components_to_samples"),
        touch(component + "/add_components_to_samples"),
    params:
        rule_name = rule_name
    run:
        rule_name = str(params.rule_name)
        run_folder = str(input.run_folder)
        component = str(input.component)

        sys.stdout.write("Started {}\n".format(rule_name))
        config = datahandling.load_config()
        unique_sample_names = {}
        for file in sorted(os.listdir(run_folder)):
            result = re.search(config["read_pattern"], file)
            if result and os.path.isfile(os.path.realpath(os.path.join(run_folder, file))):
                sample_name = result.group("sample_name")
                unique_sample_names[sample_name] = unique_sample_names.get(sample_name, 0) + 1

        for sample_name in unique_sample_names:
            sample_config = sample_name + "/sample.yaml"
            sample_db = datahandling.load_sample(sample_config)
            sample_db["components"] = sample_db.get("components", [])
            for component_name in components:
                component_id = datahandling.load_component(os.path.join(component, component_name + ".yaml")).get("_id",)
                if component_id is not None:
                    insert_component = True
                    for sample_component in sample_db["components"]:
                        if component_id == sample_component["_id"]:
                            insert_component = False
                    if insert_component is True:
                        sample_db["components"].append({"name": component_name, "_id": component_id})
            datahandling.save_sample(sample_db, sample_config)
        sys.stdout.write("Done {}\n".format(rule_name))


rule_name = "initialize_sample_components_for_each_sample"
rule initialize_sample_components_for_each_sample:
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
        rules.add_components_to_samples.output,
        component = component,
        run_folder = run_folder
    output:
        touch(rerun_folder + "/initialize_sample_components_for_each_sample"),
        touch(component + "/initialize_sample_components_for_each_sample"),
    params:
        rule_name = rule_name
    run:
        rule_name = str(params.rule_name)
        run_folder = str(input.run_folder)

        sys.stdout.write("Started {}\n".format(rule_name))
        config = datahandling.load_config()
        unique_sample_names = {}
        for file in sorted(os.listdir(run_folder)):
            result = re.search(config["read_pattern"], file)
            if result and os.path.isfile(os.path.realpath(os.path.join(run_folder, file))):
                sample_name = result.group("sample_name")
                unique_sample_names[sample_name] = unique_sample_names.get(sample_name, 0) + 1

        for sample_name in unique_sample_names:
            sample_config = sample_name + "/sample.yaml"
            sample_db = datahandling.load_sample(sample_config)

            sample_id = sample_db.get("_id",)
            component_dict = sample_db.get("components",[])
            for item in component_dict:
                component_name = item.get("name",)
                component_id = item.get("_id",)
                sample_component_path = sample_name + "/" + sample_name + "__" + component_name + ".yaml"
                sample_component_db = datahandling.load_sample_component(sample_component_path)
                sample_component_db["sample"] = {"name": sample_name, "_id": sample_id}
                sample_component_db["component"] = {"name": component_name, "_id": component_id}
                sample_component_db["status"] = "initialized"
                datahandling.save_sample_component(sample_component_db, sample_component_path)
        sys.stdout.write("Done {}\n".format(rule_name))


rule_name = "initialize_run"
rule initialize_run:
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
        rules.copy_run_info.output,
        rules.initialize_sample_components_for_each_sample.output,
        run_folder = run_folder,
        component = component
    output:
        touch(rerun_folder + "/initialize_run"),
        touch(component + "/initialize_run"),
    params:
        rule_name = rule_name
    run:
        rule_name = str(params.rule_name)
        run_folder = str(input.run_folder)
        component = str(input.component)

        sys.stdout.write("Started {}\n".format(rule_name))
        config = datahandling.load_config()
        unique_sample_names = {}

        run_db = datahandling.load_run(component + "/run.yaml")
        run_db["name"] = config.get("run_name", os.path.realpath(os.path.join(run_folder)).split("/")[-1])
        run_db["type"] = config.get("type", "default")
        for file in sorted(os.listdir(run_folder)):
            result = re.search(config["read_pattern"], file)
            if result and os.path.isfile(os.path.realpath(os.path.join(run_folder, file))):
                sample_name = result.group("sample_name")
                unique_sample_names[sample_name] = unique_sample_names.get(sample_name, 0) + 1

        run_db["samples"] = run_db.get("samples", [])
        # Todo: handle change in samples properly
        for sample_name in unique_sample_names:
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
            component_id = datahandling.load_component(os.path.join(component, component_name + ".yaml")).get("_id",)
            if component_id is not None:
                insert_component = True
                for run_components in run_db["components"]:
                    if component_id == run_components["_id"]:
                        insert_component = False
                if insert_component is True:
                    run_db["components"].append({"name": component_name, "_id": component_id})

        datahandling.save_run(run_db, component + "/run.yaml")
        sys.stdout.write("Done {}\n".format(rule_name))


rule_name = "setup_sample_components_to_run"
rule setup_sample_components_to_run:
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
        rules.initialize_run.output,
        component = component,
        run_folder = run_folder
    output:
        touch(rerun_folder + "/setup_sample_components_to_run"),
        bash_file = "run_cmd_serumqc.sh"
    params:
        rule_name = rule_name
    run:
        rule_name = str(params.rule_name)
        run_folder = str(input.run_folder)
        component = str(input.component)
        run_cmd = str(output.bash_file)
        sys.stdout.write("Started {}\n".format(rule_name))
        config = datahandling.load_config()
        unique_sample_names = {}
        for file in sorted(os.listdir(run_folder)):
            result = re.search(config["read_pattern"], file)
            if result and os.path.isfile(os.path.realpath(os.path.join(run_folder, file))):
                sample_name = result.group("sample_name")
                unique_sample_names[sample_name] = unique_sample_names.get(sample_name, 0) + 1

        with open(run_cmd, "w") as run_cmd_handle:
            for sample_name in unique_sample_names:
                current_time = datetime.datetime.now()
                with open(sample_name + "/cmd_serumqc_{}.sh".format(current_time), "w") as command:
                    command.write("#!/bin/sh\n")
                    if config["grid"] == "torque":
                        command.write("#PBS -V -d . -w . -l ncpus={},mem={}gb -N 'serumqc_{}' -W group_list={} -A {} \n".format(config["memory"], config["threads"], sample, group, group))
                    elif config["grid"] == "slurm":
                        command.write("#SBATCH --mem={}G -p {} -c {} -J 'serumqc_{}'\n".format(config["memory"], config["partition"], config["threads"], sample_name))

                    sample_config = sample_name + "/sample.yaml"
                    sample_db = datahandling.load_sample(sample_config)
                    if sample_name not in config["samples_to_ignore"] and "R1" in sample_db["reads"] and "R2" in sample_db["reads"]:
                        for component_name in components:
                            component_file = os.path.dirname(workflow.snakefile) + "/snakefiles/" + component_name + ".smk"
                            if os.path.isfile(component_file):
                                command.write("if [ -d \"{}\" ]; then rm -r {}; fi;\n".format(component_name, component_name))
                                command.write("snakemake --cores {} -s {} --config Sample={};\n".format(config["threads"], component_file, "sample.yaml"))

                                sample_component_db = datahandling.load_sample_component(sample_name + "/" + sample_name + "__" + component_name + ".yaml")
                                sample_component_db["status"] = "queued to run"
                                sample_component_db["setup_date"] = current_time
                                datahandling.save_sample_component(sample_component_db, sample_name + "/" + sample_name + "__" + component_name + ".yaml")
                            else:
                                sys.stderr.write("Error component not found:{} {}".format(component_name, component_file))
                                sample_component_db = datahandling.load_sample_component(sample_name + "/" + sample_name + "__" + component_name + ".yaml")
                                sample_component_db["status"] = "component_missing"
                                sample_component_db["setup_date"] = current_time
                                datahandling.save_sample_component(sample_component_db, sample_name + "/" + sample_name + "__" + component_name + ".yaml")

                os.chmod(os.path.join(sample_name, "cmd_serumqc_{}.sh".format(current_time)), 0o777)
                if os.path.islink(os.path.join(sample_name, "cmd_serumqc.sh")):
                    os.remove(os.path.join(sample_name, "cmd_serumqc.sh"))
                os.symlink(os.path.realpath(os.path.join(sample_name, "cmd_serumqc_{}.sh".format(current_time))), os.path.join(sample_name, "cmd_serumqc.sh"))
                run_cmd_handle.write("cd {};\n".format(sample_name))
                if config["grid"] == "torque":
                    run_cmd_handle.write("qsub cmd_serumqc.sh;\n")  # dependent on grid engine
                elif config["grid"] == "slurm":
                    run_cmd_handle.write("sbatch cmd_serumqc.sh;\n")  # dependent on grid engine
                else:
                    run_cmd_handle.write("bash cmd_serumqc.sh;\n")
                run_cmd_handle.write("cd {};\n".format(os.getcwd()))
        sys.stdout.write("Started {}\n".format(rule_name))


rule create_end_file:
    input:
        rules.setup_sample_components_to_run.output
    output:
        rules.all.input
    shell:
        """
        bash run_cmd_serumqc.sh
        touch {output}
        """
