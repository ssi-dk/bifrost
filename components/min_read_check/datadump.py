import os
from bifrostlib import datahandling

def extract_has_min_num_of_reads(datadumpObj):
    import re
    summary, results = datadumpObj.get_summary_and_results()
    file_path = os.path.join(datadumpObj.get_component_name(), "has_min_num_of_reads")
    buffer = datahandling.read_buffer(file_path)
    results[file_path] = {}
    results[file_path]["min_read_num"] = int(re.search("min_read_num:\s*([0-9]+)", buffer, re.MULTILINE).group(1))
    summary["has_min_num_of_reads"] = True
    return (summary, results)

def datadump(sample_file, component_file, sample_component_file, log):
    datadumpObj = datahandling.DatadumpSampleComponentObj(sample_file, component_file, sample_component_file, log)
    datadumpObj.retrieve_data(extract_has_min_num_of_reads)
    datadumpObj.save()


datadump(
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.sample_component_file,
    snakemake.log)
