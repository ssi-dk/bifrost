import os
from bifrostlib import datahandling


def extract_has_min_num_of_reads(SampleComponentObj):
    import re
    summary, results = SampleComponentObj.get_summary_and_results()
    file_path = os.path.join(SampleComponentObj.get_component_name(), "has_min_num_of_reads")
    buffer = datahandling.read_buffer(file_path)
    results[file_path] = {}
    results[file_path]["min_read_num"] = int(re.search("min_read_num:\s*([0-9]+)", buffer, re.MULTILINE).group(1))
    summary["has_min_num_of_reads"] = True
    return (summary, results)


def datadump(SampleComponentObj, log):
    SampleComponentObj.start_data_dump(log)
    SampleComponentObj.run_data_dump_on_function(extract_has_min_num_of_reads, log)
    SampleComponentObj.end_data_dump(log)


datadump(
    snakemake.params.SampleComponentObj,
    snakemake.log)
