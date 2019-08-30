from bifrostlib import datahandling


def extract_cge_mlst_report_and_details(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction("data.yaml")
    results[key] = datahandling.load_yaml(file_path)
    strains = []
    for mlst_db in results[key]:
        strain = results[key][mlst_db]["mlst"]["results"]["sequence_type"]
        strains.append(strain)
    results["strain"] = strains
    summary["strain"] = strains
    return (summary, results)


def generate_report(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction()
    key = sampleComponentObj.get_file_location_key("data.yaml")
    data = []
    for mlst_db in results[key]:
        strain = results[key][mlst_db]["mlst"]["results"]["sequence_type"]
        alleles = ", ".join([results[key][mlst_db]["mlst"]["results"]["allele_profile"][i]["allel: ,
        data.append({
            "db": mlst_db,
            "strain": strain,
            "alleles": alleles
        })
    return data


def datadump(sampleComponentObj, log):
    sampleComponentObj.start_data_dump(log=log)
    sampleComponentObj.run_data_dump_on_function(extract_cge_mlst_report_and_details, log=log)
    sampleComponentObj.end_data_dump(generate_report, log=log)

datadump(
    snakemake.params.sampleComponentObj,
    snakemake.log)
