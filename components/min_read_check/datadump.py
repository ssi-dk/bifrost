import os
from bifrostlib import datahandling


def extract_has_min_num_of_reads(bifrost_sample_component_object):
    import re
    summary, results = bifrost_sample_component_object.get_summary_and_results()
    file_path = os.path.join(bifrost_sample_component_object.get_component_name(), "has_min_num_of_reads")
    buffer = datahandling.read_buffer(file_path)
    results[file_path] = {}
    results[file_path]["min_read_num"] = int(re.search("min_read_num:\s*([0-9]+)", buffer, re.MULTILINE).group(1))
    summary["has_min_num_of_reads"] = True
    return (summary, results)


def datadump(bifrost_sample_component_object, log):
    bifrost_sample_component_object.start_data_dump(log)
    bifrost_sample_component_object.run_data_dump_on_function(extract_has_min_num_of_reads)
    bifrost_sample_component_object(extract_has_min_num_of_reads, log)
    bifrost_sample_component_object.end_data_dump(log)


datadump(
    snakemake.params.bifrost_sample_component_object,
    snakemake.log)
