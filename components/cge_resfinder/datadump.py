import pkg_resources
import datetime
import os
import re
import sys
import traceback
from bifrostlib import datahandling




def extract_cge_resfinder_data(db, file_path, key, temp_data):
    buffer = datahandling.load_yaml(file_path)
    db["results"][key] = buffer
    return db


def convert_summary_for_reporter(db, file_path, key, temp_data):
    resfinder_dict = db["results"]["cge_resfinder/data_resfinder_json"]["resfinder"]["results"]
    for anti_biotic_class in resfinder_dict:
        for subclass in resfinder_dict[anti_biotic_class]:
            if resfinder_dict[anti_biotic_class][subclass] != "No hit found":
                gene_dict = resfinder_dict[anti_biotic_class][subclass]
                for gene in gene_dict:
                    db["reporter"]["content"].append([gene_dict[gene]["resistance_gene"], gene_dict[gene]["coverage"], gene_dict[gene]["identity"], anti_biotic_class, gene_dict[gene]["predicted_phenotype"]])  # table rows
    return db


def script__datadump(output, sample_file, component_file, sample_component_file, log):
    try:
        output = str(output)
        log_out = str(log.out_file)
        log_err = str(log.err_file)
        sample_db = datahandling.load_sample(sample_file)
        component_db = datahandling.load_component(component_file)
        db_sample_component = datahandling.load_sample_component(sample_component_file)
        this_function_name = sys._getframe().f_code.co_name
        global GLOBAL_component_name
        GLOBAL_component_name = component_db["name"]

        datahandling.write_log(log_out, "Started {}\n".format(this_function_name))

        # Save files to DB
        datahandling.save_files_to_db(component_db["db_values_changes"]["files"], sample_component_id=db_sample_component["_id"])

        # Initialization of values, summary and reporter are also saved into the sample
        db_sample_component["summary"] = {"component": {"_id": component_db["_id"], "_date": datetime.datetime.utcnow()}}
        db_sample_component["results"] = {}
        db_sample_component["reporter"] = component_db["db_values_changes"]["sample"]["reporter"]["resistance"]

        # Data extractions
        db_sample_component = datahandling.datadump_template(extract_cge_resfinder_data, db_sample_component, file_path=os.path.join(GLOBAL_component_name, "data_resfinder.json"))
        db_sample_component = datahandling.datadump_template(convert_summary_for_reporter, db_sample_component)

        # Save to sample component
        datahandling.save_sample_component_to_file(db_sample_component, sample_component_file)
        # Save summary and reporter results into sample
        sample_db["properties"]["resistance"] = db_sample_component["summary"]
        sample_db["reporter"]["resistance"] = db_sample_component["reporter"]
        datahandling.save_sample_to_file(sample_db, sample_file)
        open(output, 'w+').close()  # touch file

    except Exception:
        datahandling.write_log(log_out, "Exception in {}\n".format(this_function_name))
        datahandling.write_log(log_err, str(traceback.format_exc()))
        raise Exception
        return 1

    finally:
        datahandling.write_log(log_out, "Done {}\n".format(this_function_name))
        return 0


script__datadump(
    snakemake.output.complete,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.sample_component_file,
    snakemake.log)
