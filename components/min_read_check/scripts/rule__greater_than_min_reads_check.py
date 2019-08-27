# script for use with snakemake
import sys
import traceback
import re
from bifrostlib import datahandling


def rule__greater_than_min_reads_check(input, output, sample_file, component_file, log):
    this_function_name = sys._getframe().f_code.co_name
    ruleObj = datahandling.snakemakeRuleScriptObj(input, output, sample_file, component_file, log, this_function_name)
    db_sample, db_component = ruleObj.get_sample_and_component_dbs()
    try:
        # Variables being used
        min_read_number = int(db_component["options"]["min_num_reads"])
        stats_data = datahandling.read_buffer(input.stats_file)
        output_file = str(output.file)

        # Code to run
        if int(re.search("Result:\s*([0-9]+)\sreads", stats_data, re.MULTILINE).group(1)) > min_read_number:
            with open(output_file, "w") as output:
                output.write("min_read_num:{}".format(min_read_number))

    except Exception:
        ruleObj.write_log_err(str(traceback.format_exc()))

    finally:
        return ruleObj.rule_done()


rule__greater_than_min_reads_check(
    snakemake.input,
    snakemake.output,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.log)
