from bifrostlib import datahandling


def extract_ariba_mlst_report_and_details(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction("data.yaml")
    results[key] = datahandling.load_yaml(file_path)
    strains = []
    for mlst_db in results[key]:
        strain_db = results[key][mlst_db]["report"]
        strain = strain_db["ST"]
        strains.append(strain)
        results["strain"] = strains
        summary["strain"] = strains
    return (summary, results)


def generate_report(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction()
    key = sampleComponentObj.get_file_location_key("data.yaml")

    data = []
    for mlst_db in results[key]:
        strain_db = results[key][mlst_db]["report"]
        alleles = []
        strain = strain_db["ST"]
        for gene in strain_db:
            if gene != "ST":
                alleles.append("{}_{}".format(gene, strain_db[gene]))
        alleles = ", ".join(alleles)
        data.append({
            "db": mlst_db,
            "strain": strain,
            "alleles": alleles
        })
    return data


def datadump(sampleComponentObj, log):
    sampleComponentObj.start_data_dump(log=log)
    sampleComponentObj.run_data_dump_on_function(extract_ariba_mlst_report_and_details, log=log)
    sampleComponentObj.end_data_dump(generate_report, log=log)

datadump(
    snakemake.params.sampleComponentObj,
    snakemake.log)
