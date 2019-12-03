# script for use with snakemake
import sys
import traceback
import os
from bifrostlib import datahandling


def rule__ariba_virulencefinder(input, output, sampleComponentObj, log):
    try:
        this_function_name = sys._getframe().f_code.co_name
        name, options, resources = sampleComponentObj.start_rule(this_function_name, log=log)

        # Variables being used
        database_path = resources["database_path"]
        reads = input.reads  # expected a tuple of read locations
        output_file = output.complete  # a file to mark success for snakemake

        # Code to run
        command = "ariba run --force {} {} {} {}".format(database_path, reads[0], reads[1], os.path.join(name, "virulence"))
        sampleComponentObj.rule_run_cmd(command, log)

        sampleComponentObj.end_rule(this_function_name, log=log)
    except Exception:
        sampleComponentObj.write_log_err(log, str(traceback.format_exc()))


rule__ariba_virulencefinder(
    snakemake.input,
    snakemake.output,
    snakemake.params.sampleComponentObj,
    snakemake.log)
