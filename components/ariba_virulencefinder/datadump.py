from bifrostlib import datahandling


def extract_ariba_virulencefinder_data(sampleComponentObj):
    import pandas
    # all the rows are turned into a subelement of one gene, all subelements are expected to be the same size
    #- Sub Function start --------------------------------------------------------------------------

    def snp_info_parser(row):
        var_info = []
        num_of_snp_elements = len(row["has_known_var"])
        for i in range(num_of_snp_elements):
            var_info.append({
                "ref_base_assembled": row["ref_base_assembled"][i],
                "pc_ident": row["pc_ident"][i],
                "ctg": row["ctg"][i],
                "ctg_len": row["ctg_len"][i],
                "ctg_cov": row["ctg_cov"][i],
                "known_var": row["known_var"][i],
                "var_type": row["var_type"][i],
                "var_seq_type": row["var_seq_type"][i],
                "known_var_change": row["known_var_change"][i],
                "has_known_var": row["has_known_var"][i],
                "ref_ctg_change": row["ref_ctg_change"][i],
                "ref_ctg_effect": row["ref_ctg_effect"][i],
                "ref_start": row["ref_start"][i],
                "ref_end": row["ref_end"][i],
                "ref_nt": row["ref_nt"][i],
                "ctg_start": row["ctg_start"][i],
                "ctg_end": row["ctg_end"][i],
                "ctg_nt": row["ctg_nt"][i],
                "smtls_total_depth": row["smtls_total_depth"][i],
                "smtls_nts": row["smtls_nts"][i],
                "smtls_nts_depth": row["smtls_nts_depth"][i]})
        return var_info
    #- Sub Function end ----------------------------------------------------------------------------
    df = pandas.read_csv(file_path, sep="\t")
    # These are all the columns which will be included in the final dataframe except for the collapsed column
    grouped_df = df.groupby(["#ariba_ref_name", "ref_name", "gene", "var_only", "flag", "reads", "cluster", "ref_len", "var_description", "free_text"])
    # Turn the columns which are not part of the final dataframe into lists, some will have 1 entry some will have multiple
    flattened_df = grouped_df.agg(list)
    # Reset the indexing for formatting purposes otherwise you layer them
    flattened_df = flattened_df.reset_index()
    # Apply the function onto the columns that are in multiples and convert data accordingly and assign that to a new variable
    flattened_df["var_info"] = flattened_df.apply(snp_info_parser, axis=1)
    # Drop the columns used to make the combined value from the function
    flattened_df = flattened_df.drop(columns=["ref_base_assembled", "pc_ident", "ctg", "ctg_len", "ctg_cov", "known_var", "var_type", "var_seq_type", "known_var_change", "has_known_var", "ref_ctg_change", "ref_ctg_effect", "ref_start", "ref_end", "ref_nt", "ctg_start", "ctg_end", "ctg_nt", "smtls_total_depth", "smtls_nts", "smtls_nts_depth"])
    # Set the reference id back to a variable for the dict to be in proper format
    flattened_df = flattened_df.set_index("#ariba_ref_name")

    results[key] = flattened_df.to_dict(orient="index")
    return (summary, results)


def generate_report(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction()
    key = sampleComponentObj.get_file_location_key("virulence/report.tsv")
    data = []
    for gene in results[key]:
        variant_count = 0
        for variant in results[key][gene]["var_info"]:
            variant_count = variant_count + 1
        data.append({
            "gene": gene,
            "coverage": round(results[key][gene].get("ref_base_assembled", 0) / results[key][gene].get("ref_len", 1), 3),
            "identity": results[key][gene]["pc_ident"],
            "variants": variant_count
        })
    return data

def datadump(sampleComponentObj, log):
    sampleComponentObj.start_data_dump(log=log)
    sampleComponentObj.run_data_dump_on_function(extract_ariba_virulencefinder_data, log=log)
    sampleComponentObj.end_data_dump(generate_report_function=generate_report, log=log)

datadump(
    snakemake.params.sampleComponentObj,
    snakemake.log)
