
#---- Templated section: start ---------------------------------------------------------------------
import sys
import subprocess
import traceback
import os
from bifrostlib import datahandling
#---- Templated section: end -----------------------------------------------------------------------
#**** Dynamic section: start ***********************************************************************
def script__run_rule_name(input, output, sample_file, component_file, folder, log):
    import pandas
#**** Dynamic section: end *************************************************************************
#---- Templated section: start ---------------------------------------------------------------------
    try:
        log_out = str(log.out_file)
        log_err = str(log.err_file)
        db_sample = datahandling.load_sample(sample_file)
        db_component = datahandling.load_component(component_file)
        this_function_name = sys._getframe().f_code.co_name

        datahandling.write_log(log_out, "Started {}\n".format(this_function_name))

#---- Templated section: end -----------------------------------------------------------------------
#**** Dynamic section: start ***********************************************************************
        command = "echo \"\{'key':'value'\}\" > {} 1> {} 2> {}".format(output.file, log_out, log_err)
        datahandling.write_log(log_out, "Running:{}".format(command))
        subprocess.Popen(command, shell=True).communicate()
#**** Dynamic section: end *************************************************************************
#---- Templated section: start ---------------------------------------------------------------------

    except Exception:
        datahandling.write_log(log_out, "Exception in {}\n".format(this_function_name))
        datahandling.write_log(log_err, str(traceback.format_exc()))

    finally:
        datahandling.write_log(log_out, "Done {}\n".format(this_function_name))
        return 0


#**** Dynamic section: start ***********************************************************************
script__run_rule_name(\
#**** Dynamic section: end *************************************************************************
    snakemake.input,
    snakemake.output,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.folder,
    snakemake.log)
#---- Templated section: end -----------------------------------------------------------------------
