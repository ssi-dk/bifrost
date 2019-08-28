# script for use with snakemake
import sys
import traceback
import re
from bifrostlib import datahandling


def rule__greater_than_min_reads_check(input, output, SampleComponentObj, log):
    try:
        this_function_name = sys._getframe().f_code.co_name
        db_sample, db_component = SampleComponentObj.start_rule(this_function_name, log)
        # Variables being used
        min_read_number = int(db_component["options"]["min_num_reads"])
        stats_data = datahandling.read_buffer(input.stats_file)
        output_file = str(output.file)

        # Code to run
        if int(re.search("Result:\s*([0-9]+)\sreads", stats_data, re.MULTILINE).group(1)) > min_read_number:
            with open(output_file, "w") as output:
                output.write("min_read_num:{}".format(min_read_number))

    except Exception:
        SampleComponentObj.write_log_err(log, str(traceback.format_exc()))

    finally:
        return SampleComponentObj.end_rule(this_function_name, log)


rule__greater_than_min_reads_check(
    snakemake.input,
    snakemake.output,
    snakemake.params.SampleComponentObj,
    snakemake.log)