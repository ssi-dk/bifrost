# script for use with snakemake
import sys
import subprocess
import traceback
from bifrostlib import datahandling


def rule__run_cge_mlst(input, output, sampleComponentObj, log):
    try:
        this_function_name = sys._getframe().f_code.co_name
        name, options, resources = sampleComponentObj.start_rule(this_function_name, log=log)

        # Variables being used
        database_path = resources["database_path"]
        reads = input.reads  # expected a tuple of read locations
        output_file = output.complete  # a file to mark success for snakemake
        species = sampleComponentObj.get_sample_properties_by_category("species_detection")["species"]

        # Code to run
        if species not in options["mlst_species_mapping"]:
            sampleComponentObj.write_log_out(log, "species {} not in mlst species\n".format(species))
            sampleComponentObj.rule_run_cmd("touch {}/no_mlst_species_DB".format(name), log)
        else:
            mlst_species = options["mlst_species_mapping"][species]
            data_dict = {}
            for mlst_entry in mlst_species:
                mlst_entry_path = "{}/{}".format(name, mlst_entry)
                sampleComponentObj.rule_run_cmd("if [ -d \"{}\" ]; then rm -r {}; fi".format(mlst_entry_path, mlst_entry_path), log)
                sampleComponentObj.rule_run_cmd("ls {} -lah".format(database_path), log)

                sampleComponentObj.rule_run_cmd("mkdir {}; mlst.py -x -matrix -s {} -p {} -mp kma -i {} {} -o {}".format(mlst_entry_path, mlst_entry, database_path, reads[0], reads[1], mlst_entry_path), log)
                data_dict[mlst_entry] = datahandling.load_yaml("{}/data.json".format(mlst_entry_path))
            datahandling.save_yaml(data_dict, output_file)

        sampleComponentObj.end_rule(this_function_name, log=log)
    except Exception:
        sampleComponentObj.write_log_err(log, str(traceback.format_exc()))


rule__run_cge_mlst(
    snakemake.input,
    snakemake.output,
    snakemake.params.sampleComponentObj,
    snakemake.log)
