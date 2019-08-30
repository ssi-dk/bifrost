# script for use with snakemake
import sys
import traceback
import os
from bifrostlib import datahandling

def rule__greater_than_min_reads_check(input, output, sampleComponentObj, log):
    import re
    try:
        this_function_name = sys._getframe().f_code.co_name
        name, options, resources = sampleComponentObj.start_rule(this_function_name, log=log)

        # Variables being used
        min_read_number = options["min_num_reads"]
        stats_data = datahandling.read_buffer(input.stats_file)
        output_file = str(output.file)

        # Code to run
        num_of_reads = int(re.search("Result:\s*([0-9]+)\sreads", stats_data, re.MULTILINE).group(1))
        if num_of_reads > min_read_number:
            with open(output_file, "w") as output:
                output.write("min_read_num:{}".format(min_read_number))
        else:
            sampleComponentObj.write_log_out(log, "Doesn't have min reads {}reads found\n".format(num_of_reads))

        sampleComponentObj.end_rule(this_function_name, log=log)
    except Exception:
        sampleComponentObj.write_log_err(log, str(traceback.format_exc()))


rule__greater_than_min_reads_check(
    snakemake.input,
    snakemake.output,
    snakemake.params.sampleComponentObj,
    snakemake.log)
