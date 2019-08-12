import pkg_resources
from ruamel.yaml import YAML
import os
import re
from bifrostlib import datahandling
import sys

config = datahandling.load_config()


def extract_bbuk_log(file_path, key, data_dict):
    buffer = datahandling.read_buffer(file_path)
    data_dict["results"][key]["input_reads_num"] = int(re.search("Input:\s*([0-9]+)\sreads", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["filtered_reads_num"] = int(re.search("Result:\s*([0-9]+)\sreads", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["input_reads_bases"] = int(re.search("Input:.*?([0-9]+)\sbases", buffer, re.MULTILINE).group(1))
    data_dict["results"][key]["filtered_reads_bases"] = int(re.search("Result:.*?([0-9]+)\sbases", buffer, re.MULTILINE).group(1))
    data_dict["summary"]["filtered_reads_num"] = data_dict["results"][key]["filtered_reads_num"]
    return data_dict


def script__datadump(folder, sample_file, component_file, sample_component_file, log):
    try:
        log_out = str(log.out_file)
        log_err = str(log.err_file)
        db_sample = datahandling.load_sample(sample_file)
        db_component = datahandling.load_component(component_file)
        db_sample_component = datahandling.load_sample_component(sample_component_file)
        this_function_name = sys._getframe().f_code.co_name

        datahandling.log(log_out, "Started {}\n".format(this_function_name))
        
        db_sample_component["summary"] = {"component": {"id": db_component["_id"], "date": datetime.datetime.utcnow()}}
        db_sample_component["results"] = {}

        db_sample_component = datahandling.datadump_template(db_sample_component, folder, "log/setup__filter_reads_with_bbduk.err.log", extract_bbuk_log)

        datahandling.save_sample_component(db_sample_component, sample_component_file)

    except Exception:
        datahandling.log(log_out, "Exception in {}\n".format(this_function_name))
        datahandling.log(log_err, str(traceback.format_exc()))

    finally:
        datahandling.log(log_out, "Done {}\n".format(this_function_name))
        return 0


script__datadump(
    snakemake.params.folder,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.sample_component_file,
    snakemake.log)
