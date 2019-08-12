# script for use with snakemake
import sys
import subprocess
import traceback
from bifrostlib import datahandling


def script__run_cge_mlst(input, output, sample_file, component_file, folder, log):
    try:
        log_out = str(log.out_file)
        log_err = str(log.err_file)
        db_sample = datahandling.load_sample(sample_file)
        db_component = datahandling.load_component(component_file)
        this_function_name = sys._getframe().f_code.co_name

        datahandling.log(log_out, "Started {}\n".format(this_function_name))

        # Variables being used
        database_path = db_component["database_path"]
        reads = input.reads  # expected a tuple of read locations
        output_file = output.complete  # a file to mark success for snakemake
        species = db_sample["properties"]["species"]

        # Code to run
        if species not in db_component["mlst_species_mapping"]:
            datahandling.log(log_out, "cge mlst species: {}\n".format(species))
            subprocess.Popen("touch " + folder + "/no_mlst_species_DB").communicate()
        else:
            mlst_species = db_component["mlst_species_mapping"][species]
            data_dict = {}
            for mlst_entry in mlst_species:
                mlst_entry_path = folder + "/" + mlst_entry
                datahandling.log(log_out, "mlst {} on species: {}\n".format(mlst_entry, species))
                subprocess.Popen("if [ -d \"{}\" ]; then rm -r {}; fi".format(mlst_entry_path, mlst_entry_path), shell=True).communicate()
                subprocess.Popen("mkdir {}; mlst.py -x -matrix -s {} -p {} -mp kma -i {} {} -o {} 1> {} 2> {}".format(mlst_entry_path, mlst_entry, database_path, reads[0], reads[1], mlst_entry_path, log_out, log_err), shell=True).communicate()
                data_dict[mlst_entry] = datahandling.load_yaml(mlst_entry + "/data.json")
            datahandling.save_yaml(data_dict, output_file)

    except Exception:
        datahandling.log(log_out, "Exception in {}\n".format(this_function_name))
        datahandling.log(log_err, str(traceback.format_exc()))

    finally:
        datahandling.log(log_out, "Done {}\n".format(this_function_name))
        return 0


script__run_cge_mlst(
    snakemake.input,
    snakemake.output,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.folder,
    snakemake.log)
