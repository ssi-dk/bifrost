
#---- Templated section: start ---------------------------------------------------------------------
import sys
import subprocess
import traceback
import os
from bifrostlib import datahandling
#---- Templated section: end -----------------------------------------------------------------------
#**** Dynamic section: start ***********************************************************************
def script__run_ariba_mlst(input, output, sample_file, component_file, folder, log):
#**** Dynamic section: end *************************************************************************
#---- Templated section: start ---------------------------------------------------------------------
    try:
        log_out = str(log.out_file)
        log_err = str(log.err_file)
        db_sample = datahandling.load_sample(sample_file)
        db_component = datahandling.load_component(component_file)
        this_function_name = sys._getframe().f_code.co_name

        datahandling.log(log_out, "Started {}\n".format(this_function_name))
#---- Templated section: end -----------------------------------------------------------------------
#**** Dynamic section: start ***********************************************************************
        # Variables being used
        database_path = db_component["database_path"]
        reads = input.reads  # expected a tuple of read locations
        output_file = output.complete  # a file to mark success for snakemake
        species = db_sample["properties"]["species"]

        # Code to run
        if species not in db_component["options"]["mlst_species_mapping"]:
            datahandling.log(log_out, "cge mlst species: {}\n".format(species))
            subprocess.Popen("touch " + folder + "/no_mlst_species_DB").communicate()
        else:
            mlst_species = db_component["options"]["mlst_species_mapping"][species]
            data_dict = {}
            for mlst_entry in mlst_species:
                mlst_entry_path = folder + "/" + mlst_entry
                mlst_database_path = os.path.join(database_path, mlst_entry)
                datahandling.log(log_out, "mlst {} on species: {}\n".format(mlst_entry, species))
                command = "if [ -d \"{}\" ]; then rm -r {}; fi".format(mlst_entry_path, mlst_entry_path)
                subprocess.Popen(command, shell=True).communicate()
                datahandling.log(log_out, "Running:{}".format(command))
                command = "mkdir {}; ariba run {} {} {} {} 1> {} 2> {}".format(mlst_entry_path, mlst_database_path, reads[0], reads[1], mlst_entry_path, log_out, log_err)
                datahandling.log(log_out, "Running:{}".format(command))
                subprocess.Popen(command, shell=True).communicate()
                data_dict[mlst_entry] = datahandling.load_yaml(folder + "/" + mlst_entry + "/data.json")
            datahandling.save_yaml(data_dict, output_file)
#**** Dynamic section: end *************************************************************************
#---- Templated section: start ---------------------------------------------------------------------
    except Exception:
        datahandling.log(log_out, "Exception in {}\n".format(this_function_name))
        datahandling.log(log_err, str(traceback.format_exc()))

    finally:
        datahandling.log(log_out, "Done {}\n".format(this_function_name))
        return 0

#**** Dynamic section: start ***********************************************************************
script__run_ariba_mlst(\
#**** Dynamic section: end *************************************************************************
    snakemake.input,
    snakemake.output,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.folder,
    snakemake.log)
#---- Templated section: end -----------------------------------------------------------------------
