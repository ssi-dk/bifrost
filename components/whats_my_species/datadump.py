import os
from bifrostlib import datahandling


def extract_bracken_txt(sampleComponentObj):
    summary, results = sampleComponentObj.get_summary_and_results()
    file_path = os.path.join(sampleComponentObj.get_component_name(), "bracken.txt")
    results[file_path] = {}
    buffer = datahandling.read_buffer(file_path)
    buffer = buffer.split("\n")
    number_of_entries = min(len(buffer) - 1, 2)
    if number_of_entries > 0:
        for i in range(1, 1 + number_of_entries):  # skip first line as it's header
            results[file_path]["species_" + str(i) + "_name"] = buffer[i].split("\t")[0]
            results[file_path]["species_" + str(i) + "_kraken_assigned_reads"] = buffer[i].split("\t")[3]
            results[file_path]["species_" + str(i) + "_added_reads"] = buffer[i].split("\t")[4]
            results[file_path]["species_" + str(i) + "_count"] = int(buffer[i].split("\t")[5].strip())
    return (summary, results)


def extract_kraken_report_bracken_txt(sampleComponentObj):
    summary, results = sampleComponentObj.get_summary_and_results()
    file_path = os.path.join(sampleComponentObj.get_component_name(), "bracken.txt")
    results[file_path] = {}
    buffer = datahandling.read_buffer(file_path)
    buffer = buffer.split("\n")
    if len(buffer) > 2:
        results[file_path]["unclassified_count"] = int(buffer[0].split("\t")[1])
        results[file_path]["root"] = int(buffer[1].split("\t")[1])
    return (summary, results)


def species_math(sampleComponentObj):
    summary, results = sampleComponentObj.get_summary_and_results()
    if "status" not in results[sampleComponentObj.get_component_name() + "/kraken_report_bracken_txt"] and "status" not in results[sampleComponentObj.get_component_name() + "/bracken_txt"] and "species_1_count" in results[sampleComponentObj.get_component_name() + "/bracken_txt"] and "species_2_count" in results[sampleComponentObj.get_component_name() + "/bracken_txt"]:
        summary["percent_unclassified"] = results[sampleComponentObj.get_component_name() + "/kraken_report_bracken_txt"]["unclassified_count"] / (results[sampleComponentObj.get_component_name() + "/kraken_report_bracken_txt"]["unclassified_count"] + results[sampleComponentObj.get_component_name() + "/kraken_report_bracken_txt"]["root"])
        summary["percent_classified_species_1"] = results[sampleComponentObj.get_component_name() + "/bracken_txt"]["species_1_count"] / (results[sampleComponentObj.get_component_name() + "/kraken_report_bracken_txt"]["unclassified_count"] + results[sampleComponentObj.get_component_name() + "/kraken_report_bracken_txt"]["root"])
        summary["name_classified_species_1"] = results[sampleComponentObj.get_component_name() + "/bracken_txt"]["species_1_name"]
        summary["percent_classified_species_2"] = results[sampleComponentObj.get_component_name() + "/bracken_txt"]["species_2_count"] / (results[sampleComponentObj.get_component_name() + "/kraken_report_bracken_txt"]["unclassified_count"] + results[sampleComponentObj.get_component_name() + "/kraken_report_bracken_txt"]["root"])
        summary["name_classified_species_2"] = results[sampleComponentObj.get_component_name() + "/bracken_txt"]["species_2_name"]
        summary["detected_species"] = summary["name_classified_species_1"]
    return (summary, results)


def set_sample_species(sampleComponentObj):
    summary, results = sampleComponentObj.get_summary_and_results()
    sample_info = sampleComponentObj.get_sample_properties_by_category("sample_info")
    if sample_info["provided_species"] is not None:
        summary["species"] = sample_info["provided_species"]
    else:
        summary["species"] = summary["detected_species"]
    return (summary, results)


def datadump(sampleComponentObj, log):
    sampleComponentObj.start_data_dump(log=log)
    sampleComponentObj.run_data_dump_on_function(extract_bracken_txt, log=log)
    sampleComponentObj.run_data_dump_on_function(extract_kraken_report_bracken_txt, log=log)
    sampleComponentObj.run_data_dump_on_function(species_math, log=log)
    sampleComponentObj.end_data_dump(log=log)


datadump(
    snakemake.params.sampleComponentObj,
    snakemake.log)
