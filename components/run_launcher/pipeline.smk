#!/usr/bin/env python3
"""
Initialization program for paired end Illumina reads
"""
import re
import os
import pandas
import traceback
from bifrostlib import datahandling

os.umask(0o2)

# Saving the config

datahandling.save_yaml(config, "config.yaml")
components = config["components"].split(",")
# raw_data_folder = config["raw_data_folder"]
# rename_samples = config["rename_samples"]
sample_folder = config["sample_folder"]
sample_sheet = config["sample_sheet"]


# Find samples
def find_samples_in_dir(input_folder: str, regex_pattern: str, sample_sheet: str):
options["samples_to_include"]
options["samples_to_exclude"]
options["run_directory"]
options["sample_sheet"]

samples = []


os.chdir("/Users/kimn/playground/fake_read_names")
input_folder = "."
all_items_in_dir = os.listdir(input_folder)
regex_pattern = "^(?P<sample_name>[a-zA-Z0-9\_\-]+?)(_S[0-9]+)?(_L[0-9]+)?_(R?)(?P<paired_read_number>[1|2])(_[0-9]+)?(\.fastq\.gz)$"
samples = [(i, re.search(regex_pattern,i).group("sample_name"),  re.search(regex_pattern,i).group("paired_read_number")) for i in all_items_in_dir if re.search(regex_pattern,i)]
samples.sort()
in_dir = set(all_items_in_dir)
sample_dict = {}

for item in samples:
    in_dir.remove(item[0])

for item in samples:
    sample_dict[item[1]] = []

for item in samples:
    sample_dict[item[1]].append(item[0])

unused_files = []
for item in samples:
    if len(sample_dict[item[1]]) != 2:
        in_dir.add(item[0])
        sample_dict.pop(item[1])

unused_files = list(in_dir)
unused_files.sort()

sample_sheet = "run_metadata.txt"

if os.path.isfile(sample_sheet):
    if sample_sheet in unused_files:
        unused_files.pop(unused_files.index(sample_sheet))

df = pandas.read_table(sample_sheet)
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
valid_sample_ids = list(set(df[sample_key].to_list()))
df["haveReads"] = False

for sample in sample_dict:
    if sample in valid_sample_ids:
        df.at[df[df[sample_key] == sample].index[0], "haveReads"] = True

samples = []
for sample in sample_dict:
    if sample in valid_sample_ids:
        sampleObj = datahandling.Sample(name=sample)
        datafiles = datahandling.Category(name="datafiles")
        datafiles.set_summary({"paired_data": [sample_dict[sample][0], sample_dict[sample][1]]})
        sampleObj.set_properties_datafiles(datafiles)
        sample_info = datahandling.Category(name="sample_info")
        sample_info.set_summary(df.iloc[df[df[sample_key] == sample].index[0]].to_dict())
        sampleObj.set_properties_sample_info(sample_info)
        # sample.save()
        print(sampleObj.display())
        samples.append(sampleObj)

run = datahandling.Run(name=os.getcwd().split("/")[-1])
run.set_type = "routine"
run.set_path = os.getcwd()
run.set_samples(samples)
run.set_issues(
    duplicate_samples = list(df[df['duplicatedSampleIDs']==True]['SampleID']),
    modified_samples = list(df[df['changedSampleIDs']==True]['SampleID']),
    unused_files = unused_files,
    samples_without_reads = list(df[df['haveReads']==True]['SampleID'])
)
run.set_comments("Hello")
# run.save()
print(run.display())
# run.save()

pre_script = """;
echo "start pre_script $run.name $run.type";

echo "end pre_script";

"""
positions_to_replace = re.findall(re.compile("\$run.[a-zA-Z]+"), pre_script)
for item in positions_to_replace:
    (key, value) = (item.split("."))
    pre_script = pre_script.replace(item, run.get(value))

print(pre_script)

per_sample_script = """
# Run against each sample
echo "Running $sample.name from $run.name"
singularity run --userns -B \
$sample.properties.datafiles.summary.paired_data[0],$sample.properties.datafiles.summary.paired_data[1] \
/home/projects/ssi_10003/apps/singularity/bifrost-ssi_stamper_2.0.5 \
-id $sample._id
"""
for sample in samples:
    positions_to_replace = re.findall(re.compile("\$run.[\._a-zA-Z]+"), per_sample_script)
    for item in positions_to_replace:
        (key, value) = (item.split("."))
        per_sample_script = per_sample_script.replace(item, run.get(value))
    positions_to_replace = re.findall(re.compile("\$sample\.[\.\[\]_a-zA-Z0-9]+"), per_sample_script)
    for item in positions_to_replace:
        (item.split(".")[1:])
        level = sample.display()
        for value in item.split(".")[1:]:
            if value.endswith("]"):
                (array_item, index) = value.split("[")
                index = int(index[:-1])
                print(array_item, index)
                level = level[array_item][index]
            else:
                level = level[value]
        
        per_sample_script = per_sample_script.replace(item, level)
        # (key, value) = (item.split("."))
        # print(value)
        # pre_script = pre_script.replace(item, sample.get(value))

post_script = """
"""

for sample in sample_dict:
    
    datahandling.SampleObj()

    print(sample_dict)

    print(in_dir)

#run pre script
# any special handling should be done via environmental variables
#run script
#run post script

