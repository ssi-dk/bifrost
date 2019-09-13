from bifrostlib import datahandling


def extract_has_min_num_of_reads(sampleComponentObj):
    import re
    summary, results, file_path, key = sampleComponentObj.start_data_extraction("has_min_num_of_reads")
    buffer = datahandling.read_buffer(file_path)
    results[key]["has_min_num_of_reads"] = "has_min_num_of_reads:TRUE" in buffer
    results[key]["min_read_num"] = int(re.search("min_read_num:\s*([0-9]+)", buffer, re.MULTILINE).group(1))
    summary["has_min_num_of_reads"] = results[key]["has_min_num_of_reads"]
    return (summary, results)


def datadump(sampleComponentObj, log):
    sampleComponentObj.start_data_dump(log=log)
    sampleComponentObj.run_data_dump_on_function(extract_has_min_num_of_reads, log=log)
    sampleComponentObj.end_data_dump(log=log)


datadump(
    snakemake.params.sampleComponentObj,
    snakemake.log)
