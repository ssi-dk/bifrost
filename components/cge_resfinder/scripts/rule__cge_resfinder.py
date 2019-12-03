# script for use with snakemake
import sys
import subprocess
import traceback
from bifrostlib import datahandling


def rule__run_cge_resfinder(input, output, sampleComponentObj, log):
    try:
        this_function_name = sys._getframe().f_code.co_name
        name, options, resources = sampleComponentObj.start_rule(this_function_name, log=log)

        # Variables being used
        database_path = resources["database_path"]
        reads = input.reads  # expected a tuple of read locations

        # Code to run
        sampleComponentObj.rule_run_cmd("resfinder.py -x -matrix -p {} -mp kma -i {} {} -o {}".format(
            database_path, reads[0], reads[1], name), log)

        sampleComponentObj.end_rule(this_function_name, log=log)
    except Exception:
        sampleComponentObj.write_log_err(log, str(traceback.format_exc()))


rule__run_cge_resfinder(
    snakemake.input,
    snakemake.output,
    snakemake.params.sampleComponentObj,
    snakemake.log)
