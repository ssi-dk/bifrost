from bifrostlib import datahandling


def extract_bracken_txt(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction("bracken.txt")
    buffer = datahandling.read_buffer(file_path)
    buffer = buffer.split("\n")
    number_of_entries = min(len(buffer) - 1, 2)
    if number_of_entries > 0:
        for i in range(1, 1 + number_of_entries):  # skip first line as it's header
            results[key]["species_" + str(i) + "_name"] = buffer[i].split("\t")[0]
            results[key]["species_" + str(i) + "_kraken_assigned_reads"] = buffer[i].split("\t")[3]
            results[key]["species_" + str(i) + "_added_reads"] = buffer[i].split("\t")[4]
            results[key]["species_" + str(i) + "_count"] = int(buffer[i].split("\t")[5].strip())
    return (summary, results)


def extract_kraken_report_bracken_txt(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction("kraken_report_bracken.txt")
    buffer = datahandling.read_buffer(file_path)
    buffer = buffer.split("\n")
    if len(buffer) > 2:
        results[key]["unclassified_count"] = int(buffer[0].split("\t")[1])
        results[key]["root"] = int(buffer[1].split("\t")[1])
    return (summary, results)


def species_math(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction()
    kraken_report_bracken_key = sampleComponentObj.get_file_location_key("kraken_report_bracken.txt")
    bracken_key = sampleComponentObj.get_file_location_key("bracken.txt")
    if ("status" not in results[kraken_report_bracken_key] and 
        "status" not in results[bracken_key] and 
        "species_1_count" in results[bracken_key] and 
        "species_2_count" in results[bracken_key]):
        summary["percent_unclassified"] = results[kraken_report_bracken_key]["unclassified_count"] / (results[kraken_report_bracken_key]["unclassified_count"] + results[kraken_report_bracken_key]["root"])
        summary["percent_classified_species_1"] = results[bracken_key]["species_1_count"] / (results[kraken_report_bracken_key]["unclassified_count"] + results[kraken_report_bracken_key]["root"])
        summary["name_classified_species_1"] = results[bracken_key]["species_1_name"]
        summary["percent_classified_species_2"] = results[bracken_key]["species_2_count"] / (results[kraken_report_bracken_key]["unclassified_count"] + results[kraken_report_bracken_key]["root"])
        summary["name_classified_species_2"] = results[bracken_key]["species_2_name"]
        summary["detected_species"] = summary["name_classified_species_1"]
    return (summary, results)


def set_sample_species(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction()
    sample_info = sampleComponentObj.get_sample_properties_by_category("sample_info")
    if sample_info.get("provided_species", None) is not None:
        summary["species"] = sample_info["provided_species"]
    else:
        summary["species"] = summary.get("detected_species", None)
    return (summary, results)


def datadump(sampleComponentObj, log):
    sampleComponentObj.start_data_dump(log=log)
    sampleComponentObj.run_data_dump_on_function(extract_bracken_txt, log=log)
    sampleComponentObj.run_data_dump_on_function(extract_kraken_report_bracken_txt, log=log)
    sampleComponentObj.run_data_dump_on_function(species_math, log=log)
    sampleComponentObj.run_data_dump_on_function(set_sample_species, log=log)
    sampleComponentObj.end_data_dump(log=log)


datadump(
    snakemake.params.sampleComponentObj,
    snakemake.log)
