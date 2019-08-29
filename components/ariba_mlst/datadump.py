import os
from bifrostlib import datahandling


def extract_mlst_report_and_details(sampleComponentObj):
    summary, results = sampleComponentObj.get_summary_and_results()
    file_path = os.path.join(sampleComponentObj.get_component_name(), "data.yaml")
    buffer = datahandling.load_yaml(file_path)
    key = file_path.replace(".", "_").replace("$", ".")
    results[key] = buffer
    strains = []
    for mlst_db in results[key]:
        strain_db = results[key][mlst_db]["report"]
        strain = strain_db["ST"]
        strains.append(strain)
        results["strain"] = strains
        summary["strain"] = strains
    return (summary, results)


def convert_summary_for_reporter(sampleComponentObj):
    summary, results = sampleComponentObj.get_summary_and_results()
    strains = []
    data = []
    component_name = sampleComponentObj.get_component_name()
    for mlst_db in results[component_name + "/data_yaml"]:
        strain_db = results[component_name +
                            "/data_yaml"][mlst_db]["report"]
        alleles = []
        strain = strain_db["ST"]
        strains.append(strain)
        for gene in strain_db:
            if gene != "ST":
                alleles.append("{}_{}".format(gene, strain_db[gene]))
        alleles = ", ".join(alleles)
        data.append([mlst_db, strain, alleles])
    return data


def datadump(sampleComponentObj, log):
    sampleComponentObj.start_data_dump(log=log)
    sampleComponentObj.run_data_dump_on_function(extract_mlst_report_and_details, log=log)
    sampleComponentObj.end_data_dump(convert_summary_for_reporter, log=log)

datadump(
    snakemake.params.sampleComponentObj,
    snakemake.log)
