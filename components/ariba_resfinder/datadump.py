import pkg_resources
import datetime
import os
import re
import sys
import traceback
from bifrostlib import datahandling




def test(row):
    l = len(row["ctg_start"])
    out = []
    for i in range(l):
        out.append({'start': row["ref_start"][i], 'end': row["ref_end"][i]})
    return out

def extract_ariba_resfinder_data(db, file_path, key, temp_data):
    import pandas

    # all the rows are turned into a subelement of one gene, all subelements are expected to be the same size
    def snp_info_parser(row):
        var_info = []
        num_of_snp_elements = len(row["has_known_var"])
        for i in range(num_of_snp_elements):
            var_info.append({
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

    df = pandas.read_csv(file_path, sep="\t")
    # These are all the columns which will be included in the final dataframe except for the collapsed column
    grouped_df = df.groupby(["#ariba_ref_name", "ref_name", "gene", "var_only", "flag", "reads", "cluster", "ref_len", "ref_base_assembled", "pc_ident", "ctg", "ctg_len", "ctg_cov", "var_description", "free_text"])
    # Turn the columns which are not part of the final dataframe into lists, some will have 1 entry some will have multiple
    flattened_df = grouped_df.agg(lambda x: list(x))
    # Reset the indexing for formatting purposes otherwise you layer them
    flattened_df = flattened_df.reset_index()
    # Apply the function onto the columns that are in multiples and convert data accordingly and assign that to a new variable
    flattened_df["var_info"] = flattened_df.apply(snp_info_parser, axis=1)
    # Drop the columns used to make the combined value from the function
    flattened_df = flattened_df.drop(columns=["known_var", "var_type", "var_seq_type", "known_var_change", "has_known_var", "ref_ctg_change", "ref_ctg_effect", "ref_start", "ref_end", "ref_nt", "ctg_start", "ctg_end", "ctg_nt", "smtls_total_depth", "smtls_nts", "smtls_nts_depth"])
    # Set the reference id back to a variable for the dict to be in proper format
    flattened_df = flattened_df.set_index("#ariba_ref_name")

    db["results"][key] = flattened_df.to_dict(orient="index")
    return db


def convert_summary_for_reporter(db, file_path, key, temp_data):
    report_results = db["results"]["ariba_resfinder/resistance/report_tsv"]
    for gene in report_results:
        variant_count = 0
        for variant in report_results[gene]["var_info"]:
            variant_count = variant_count + 1
        db["reporter"]["data"].append({
            "gene": gene,
            "coverage": round(report_results[gene].get("ref_base_assembled", 0) / report_results[gene].get("ref_len", 1), 3),
            "identity": report_results[gene]["pc_ident"],
            "variants": variant_count
        })
    return db


def script__datadump(output, sample_file, component_file, sample_component_file, log):
    try:
        output = str(output)
        log_out = str(log.out_file)
        log_err = str(log.err_file)
        db_sample = datahandling.load_sample(sample_file)
        db_component = datahandling.load_component(component_file)
        db_sample_component = datahandling.load_sample_component(sample_component_file)
        this_function_name = sys._getframe().f_code.co_name
        global GLOBAL_component_name
        GLOBAL_component_name = db_component["name"]
        global GLOBAL_category_name
        GLOBAL_category_name = db_component["category"]


        datahandling.log(log_out, "Started {}\n".format(this_function_name))

        # Save files to DB
        datahandling.save_files_to_db(db_component["db_values_changes"]["files"], sample_component_id=db_sample_component["_id"])

        # Initialization of values, summary and reporter are also saved into the sample
        db_sample_component["summary"] = {"component": {"_id": db_component["_id"], "_date": datetime.datetime.utcnow()}}
        db_sample_component["results"] = {}
        db_sample_component["reporter"] = db_component["db_values_changes"]["sample"]["reporter"][GLOBAL_category_name]

        # Data extractions
        db_sample_component = datahandling.datadump_template(extract_ariba_resfinder_data, db_sample_component, file_path=os.path.join(GLOBAL_component_name, "resistance/report.tsv"))
        db_sample_component = datahandling.datadump_template(convert_summary_for_reporter, db_sample_component)

        # Save to sample component
        datahandling.save_sample_component(db_sample_component, sample_component_file)
        # Save summary and reporter results into sample
        db_sample["reporter"][GLOBAL_category_name] = db_sample_component["reporter"]
        datahandling.save_sample(db_sample, sample_file)

        open(output, 'w+').close()  # touch file

    except Exception:
        datahandling.log(log_out, "Exception in {}\n".format(this_function_name))
        datahandling.log(log_err, str(traceback.format_exc()))
        raise Exception
        return 1

    finally:
        datahandling.log(log_out, "Done {}\n".format(this_function_name))
        return 0


script__datadump(
    snakemake.output.complete,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.sample_component_file,
    snakemake.log)
