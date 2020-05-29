#!/usr/bin/env python3
"""
Initialization program for paired end Illumina reads
"""
import re
import sys
import os
import numpy
import pandas
import traceback
from bifrostlib import datahandling
import pprint

os.umask(0o2)
pp = pprint.PrettyPrinter(indent=4)


def initialize_run(input_folder: str = ".", run_metadata: str = "run_metadata.txt", regex_pattern: str ="^(?P<sample_name>[a-zA-Z0-9\_\-]+?)(_S[0-9]+)?(_L[0-9]+)?_(R?)(?P<paired_read_number>[1|2])(_[0-9]+)?(\.fastq\.gz)$") -> object:
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
    sample_key = "SampleID"
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    samples_no_index = df[df[sample_key].isna()].index
    df = df.drop(samples_no_index)
    df[sample_key] = df[sample_key].astype('str')
    df["tempSampleID"] = df[sample_key]
    df[sample_key] = df[sample_key].apply(lambda x: x.strip())
    df[sample_key] = df[sample_key].str.replace(re.compile("[^a-zA-Z0-9\-\_]"),"_")
    df["changedSampleIDs"] = df['SampleID'] != df['tempSampleID']
    df["duplicatedSampleIDs"] = df.duplicated(subset=sample_key,keep="first")
    valid_sample_names = list(set(df[sample_key].tolist()))
    df["haveReads"] = False
    df["haveMetaData"] = True

    samples = []
    run = datahandling.Run(name=os.getcwd().split("/")[-1])

    for sample in sample_dict:
        if sample in valid_sample_names:
            df.loc[df[sample_key] == sample, "haveMetaData"] = True
            df.loc[df[sample_key] == sample, "haveReads"] = True
            sampleObj = datahandling.Sample(name=sample)
            datafiles = datahandling.Category(name="datafiles")
            datafiles.set_summary({"paired_data": [sample_dict[sample][0], sample_dict[sample][1]]})
            sampleObj.set_properties_datafiles(datafiles)
            sample_info = datahandling.Category(name="sample_info")
            metadata_dict = df.iloc[df[df[sample_key] == sample].index[0]].to_dict()
            for key in metadata_dict:
                if type(metadata_dict[key]) == numpy.bool_:
                    metadata_dict[key] = bool(metadata_dict[key])
            sample_info.set_summary(df.iloc[df[df[sample_key] == sample].index[0]].to_dict())
            sampleObj.set_properties_sample_info(sample_info)
            sampleObj.save()
            # pp.pprint(sampleObj.display())
            samples.append(sampleObj)
        else:
            new_row_df = pandas.DataFrame({'SampleID':[sample], 'haveReads':[True], 'haveMetaData':[False]})
            df = df.append(new_row_df, ignore_index=True, sort=False)

    run = datahandling.Run(name=os.getcwd().split("/")[-1])
    run.set_type = "routine"
    run.set_path = os.getcwd()
    run.set_samples(samples)
    run.set_issues(
        duplicate_samples = list(df[df['duplicatedSampleIDs']==True]['SampleID']),
        modified_samples = list(df[df['changedSampleIDs']==True]['SampleID']),
        unused_files = unused_files,
        samples_without_reads = list(df[df['haveReads']==True]['SampleID']),
        samples_without_metadata = list(df[df['haveMetaData']==False]['SampleID'])
    )
    run.set_comments("Hello")
    # Note when you save the run you create the ID's
    run.save() 
    pp.pprint(run.display())
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


def main(argv) -> None:
    run, samples = initialize_run(input_folder = argv[1], run_metadata = argv[2])
    script =generate_run_script(run, samples, argv[3], argv[4], argv[5])
    print(script)


if __name__ == "__main__":
    main(sys.argv)
