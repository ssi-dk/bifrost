import pkg_resources
import datetime
import os
import re
import sys
import traceback
from bifrostlib import datahandling




def extract_has_min_num_of_reads(db, file_path, key, temp_data):
    buffer = datahandling.read_buffer(file_path)
    db["results"][key]["min_read_num"] = int(re.search("min_read_num:\s*([0-9]+)", buffer, re.MULTILINE).group(1))
    db["summary"]["has_min_num_of_reads"] = True
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

        datahandling.log(log_out, "Started {}\n".format(this_function_name))

        # Save files to DB
        datahandling.save_files_to_db(db_component["db_values_changes"]["files"], sample_component_id=db_sample_component["_id"])

        # Initialization of values, summary and reporter are also saved into the sample
        db_sample_component["summary"] = {"component": {"_id": db_component["_id"], "_date": datetime.datetime.utcnow()}}
        db_sample_component["results"] = {}
        db_sample_component["reporter"] = {} # Currently unused, set to dict of component config path when used

        # Data extractions
        db_sample_component = datahandling.datadump_template(extract_has_min_num_of_reads, db_sample_component, file_path=os.path.join(GLOBAL_component_name, "has_min_num_of_reads"))

        # Save to sample component
        datahandling.save_sample_component(db_sample_component, sample_component_file)
        # Save summary and reporter results into sample
        db_sample["properties"]["size_check"] = db_sample_component["summary"]
        # db_sample["reporter"]["size_check"] = db_sample_component["reporter"]
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
