import os
from bifrostlib import datahandling


def extract_has_min_num_of_reads(sampleComponentObj):
    import re
    summary, results = sampleComponentObj.get_summary_and_results()
    file_path = os.path.join(sampleComponentObj.get_component_name(), "has_min_num_of_reads")
    key = file_path.replace(".","_").replace("$","_")
    
    buffer = datahandling.read_buffer(file_path)
    results[key] = {}
    results["FAKEKEY"]["min_read_num"] = int(re.search("min_read_num:\s*([0-9]+)", buffer, re.MULTILINE).group(1))
    summary["has_min_num_of_reads"] = True
    return (summary, results)


def datadump(sampleComponentObj, log):
    sampleComponentObj.start_data_dump(log=log)
    sampleComponentObj.run_data_dump_on_function(extract_has_min_num_of_reads, log=log)
    sampleComponentObj.end_data_dump(log=log)


datadump(
    snakemake.params.sampleComponentObj,
    snakemake.log)
