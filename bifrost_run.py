#!/usr/bin/env python3
"""
Initialization program for paired end Illumina reads
"""
import argparse
import re
import sys
import os
import numpy
import pandas
import traceback
import json
import subprocess
from bifrostlib import datahandling
from bifrostlib import mongo_interface
import pprint

os.umask(0o2)
pp = pprint.PrettyPrinter(indent=4)


def parse_args() -> object:
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument('-pre', '--pre_script',
                        help='Pre script template run before sample script')
    parser.add_argument('-per', '--per_sample_script',
                        help='Per sample script template run on each sample')
    parser.add_argument('-post', '--post_script',
                        help='Post script template run after sample script')
    parser.add_argument('-meta', '--run_metadata',
                        required=True,
                        help='Run metadata tsv')
    parser.add_argument('-reads', '--reads_folder',
                        required=True,
                        help='Run metadata tsv')
    parser.add_argument('-name', '--run_name',
                        default=None,
                        help='Run name, if not provided it will default to current folder name')
    parser.add_argument('-type', '--run_type',
                        default=None,
                        help='Run type for metadata organization')
    parser.add_argument('-metamap', '--run_metadata_column_remap',
                        help='Remaps metadata tsv columns to bifrost values')
    args: argparse.Namespace = parser.parse_args()

    setup_run(args)


def initialize_run(run_name: str, input_folder: str = ".", run_metadata: str = "run_metadata.txt", run_type: str = None, rename_column_file = None, regex_pattern: str ="^(?P<sample_name>[a-zA-Z0-9\_\-]+?)(_S[0-9]+)?(_L[0-9]+)?_(R?)(?P<paired_read_number>[1|2])(_[0-9]+)?(\.fastq\.gz)$") -> object:
    all_items_in_dir = os.listdir(input_folder)
    potential_samples = [(i, re.search(regex_pattern,i).group("sample_name"),  re.search(regex_pattern,i).group("paired_read_number")) for i in all_items_in_dir if re.search(regex_pattern,i)]
    potential_samples.sort()
    in_dir = set(all_items_in_dir)
    sample_dict = {}

    for item in potential_samples:
        in_dir.remove(item[0])

    for item in potential_samples:
        sample_dict[item[1]] = []

    for item in potential_samples:
        sample_dict[item[1]].append(item[0])

    unused_files = []
    for item in potential_samples:
        if len(sample_dict[item[1]]) != 2:
            in_dir.add(item[0])
            sample_dict.pop(item[1])

    unused_files = list(in_dir)
    unused_files.sort()

    if os.path.isfile(run_metadata):
        if run_metadata in unused_files:
            unused_files.pop(unused_files.index(run_metadata))

    df = pandas.read_table(run_metadata)
    if rename_column_file != None:
        with open(rename_column_file, "r") as rename_file:
            df = df.rename(columns=json.load(rename_file))
    sample_key = "sample_name"
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    samples_no_index = df[df[sample_key].isna()].index
    df = df.drop(samples_no_index)
    df[sample_key] = df[sample_key].astype('str')
    df["temp_sample_name"] = df[sample_key]
    df[sample_key] = df[sample_key].apply(lambda x: x.strip())
    df[sample_key] = df[sample_key].str.replace(re.compile("[^a-zA-Z0-9\-\_]"),"_")
    df["changed_sample_names"] = df['sample_name'] != df['temp_sample_name']
    df["duplicated_sample_names"] = df.duplicated(subset=sample_key,keep="first")
    valid_sample_names = list(set(df[sample_key].tolist()))
    df["haveReads"] = False
    df["haveMetaData"] = True

    samples = []
    run = datahandling.Run(name=run_name)

    for sample in sample_dict:
        if sample in valid_sample_names:
            df.loc[df[sample_key] == sample, "haveMetaData"] = True
            df.loc[df[sample_key] == sample, "haveReads"] = True
            sampleObj = datahandling.Sample(name=sample)
            datafiles = datahandling.Category(name="paired_reads")
            datafiles.set_summary({"data": [os.path.abspath(os.path.join(input_folder, sample_dict[sample][0])), os.path.abspath(os.path.join(input_folder, sample_dict[sample][1]))]})
            sampleObj.set_properties_paired_reads(datafiles)
            sample_info = datahandling.Category(name="sample_info")
            # This statement is a bit hacky, it converts to json and decodes json. This is chosen over .to_dict as that maintains numpy datatypes and I want python data types
            metadata_dict = json.loads(df.iloc[df[df[sample_key] == sample].index[0]].to_json())
            sample_info.set_summary(metadata_dict)
            sampleObj.set_properties_sample_info(sample_info)
            sampleObj.save()
            # pp.pprint(sampleObj.display())
            samples.append(sampleObj)
        else:
            new_row_df = pandas.DataFrame({'sample_name':[sample], 'haveReads':[True], 'haveMetaData':[False]})
            df = df.append(new_row_df, ignore_index=True, sort=False)

    run.set_type(run_type)
    run.set_path(os.getcwd())
    run.set_samples(samples)
    run.set_issues(
        duplicate_samples = list(df[df['duplicated_sample_names']==True]['sample_name']),
        modified_samples = list(df[df['changed_sample_names']==True]['sample_name']),
        unused_files = unused_files,
        samples_without_reads = list(df[df['haveReads']==True]['sample_name']),
        samples_without_metadata = list(df[df['haveMetaData']==False]['sample_name'])
    )
    # Note when you save the run you create the ID's
    run.save() 
    # pp.pprint(run.display())
    # df.to_csv("test.txt")

    return (run, samples)


def replace_run_info_in_script(script: str, run: object) -> str:
    positions_to_replace = re.findall(re.compile("\$run.[a-zA-Z]+"), script)
    for item in positions_to_replace:
        (key, value) = (item.split("."))
        script = script.replace(item, run.get(value))
    return script


def replace_sample_info_in_script(script: str, sample: object) -> str:
    positions_to_replace = re.findall(re.compile("\$sample\.[\.\[\]_a-zA-Z0-9]+"), script)
    for item in positions_to_replace:
        (item.split(".")[1:])
        level = sample.display()
        for value in item.split(".")[1:]:
            if value.endswith("]"):
                (array_item, index) = value.split("[")
                index = int(index[:-1])
                level = level[array_item][index]
            else:
                level = level[value]
        if(level is not None):
            if not isinstance(level, str):
                level = str(level)
            script = script.replace(item, level)
    return script


def generate_run_script(run: object, samples: object, pre_script_location: str, per_sample_script_location: str, post_script_location: str) -> str:
    script = ""
    if pre_script_location != None:
        with open(pre_script_location, "r") as pre_script_file:
            pre_script = pre_script_file.read()
            script = script + replace_run_info_in_script(pre_script, run)

    if per_sample_script_location != None:
        with open(per_sample_script_location, "r") as per_sample_script_file:
            per_sample_script = per_sample_script_file.read()
        per_sample_script = replace_run_info_in_script(per_sample_script, run)
        for sample in samples:
            script = script + replace_sample_info_in_script(per_sample_script, sample)

    if post_script_location != None:
        with open(post_script_location, "r") as post_script_file:
            post_script = post_script_file.read()
        script = script + replace_run_info_in_script(post_script, run)

    return script


def setup_run(args: object) -> str:
    run_name = args.run_name
    if run_name is None:
        run_name = os.getcwd().split("/")[-1]
    runs = mongo_interface.get_runs(names=[run_name])
    print(run_name)
    if len(runs) > 0:
        print(run_name+" already in DB, please correct before attempting to run again")
    else:
        run, samples = initialize_run(
            run_name,
            input_folder=args.reads_folder,
            run_metadata=args.run_metadata,
            rename_column_file=args.run_metadata_column_remap,
            run_type=args.run_type)
        script = generate_run_script(
            run,
            samples,
            args.pre_script,
            args.per_sample_script,
            args.post_script)
        with open("run_script.sh", "w") as fh:
            fh.write(script)
        with open("run.yaml", "w") as fh:
            fh.write(pprint.pformat(run.display()))
        with open("samples.yaml", "w") as fh:
            for sample in samples:
                fh.write(pprint.pformat(sample.display()))


if __name__ == "__main__":
    parse_args()
