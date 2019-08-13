# script for use with snakemake
import sys
import subprocess
import traceback
from bifrostlib import datahandling


def script__run_cge_resfinder(input, output, sample_file, component_file, folder, log):
    try:
        log_out = str(log.out_file)
        log_err = str(log.err_file)
        db_sample = datahandling.load_sample(sample_file)
        db_component = datahandling.load_component(component_file)
        this_function_name = sys._getframe().f_code.co_name

        datahandling.log(log_out, "Started {}\n".format(this_function_name))

        # Variables being used
        min_read_number = int(db_component["min_num_reads"])
        stats_data = datahandling.read_buffer(input.stats_file)
        output_file = str(output.file)

        # Code to run
        if int(re.search("Result:\s*([0-9]+)\sreads", stats_data, re.MULTILINE).group(1)) > min_read_number:
            with open(output_file, "w") as output:
                output.write("min_read_num:{}".format(min_read_number))

    except Exception:
        datahandling.log(log_out, "Exception in {}\n".format(this_function_name))
        datahandling.log(log_err, str(traceback.format_exc()))

    finally:
        datahandling.log(log_out, "Done {}\n".format(this_function_name))
        return 0


script__run_cge_resfinder(
    snakemake.input,
    snakemake.output,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.folder,
    snakemake.log)
