# script for use with snakemake
import sys
import traceback
import os
from bifrostlib import datahandling

def rule__run_ariba_mlst(input, output, sampleComponentObj, log):
    import pandas
    try:
        this_function_name = sys._getframe().f_code.co_name
        sampleComponentObj.start_rule(this_function_name, log=log)

        # Variables being used
        reads = input.reads  # expected a tuple of read locations
        output_file = str(output.complete)  # a file to mark success for snakemake
        database_path = sampleComponentObj.get_resources()["database_path"]
        species = sampleComponentObj.get_sample_properties_by_category("species_detection")["summary"]["species"]
        mlst_species_mapping = sampleComponentObj.get_options()["mlst_species_mapping"]

        # Code to run
        if species not in mlst_species_mapping:
            sampleComponentObj.write_log_out(log, "species {} not in mlst species\n".format(species))
            sampleComponentObj.rule_run_cmd("touch " + folder + "/no_mlst_species_DB", log)
        else:
            mlst_species = mlst_species_mapping[species]
            data_dict = {}
            for mlst_entry in mlst_species:
                data_dict[mlst_entry] = {}
                mlst_entry_path = folder + "/" + mlst_entry
                mlst_database_path = os.path.join(database_path, mlst_entry, "ref_db")
                sampleComponentObj.rule_run_cmd("ariba run --force {} {} {} {} 1> {} 2> {}".format(mlst_database_path, reads[0], reads[1], mlst_entry_path, log_out, log_err), log)
                data_dict[mlst_entry]["report"] = pandas.read_csv(os.path.join(mlst_entry_path, "mlst_report.tsv"), sep="\t").to_dict(orient="records")[0]
                data_dict[mlst_entry]["report_details"] = pandas.read_csv(os.path.join(mlst_entry_path, "mlst_report.details.tsv"), sep="\t", index_col="gene").to_dict(orient="index")
            datahandling.save_yaml(data_dict, output_file)

        sampleComponentObj.end_rule(this_function_name, log=log)
    except Exception:
        sampleComponentObj.write_log_err(log, str(traceback.format_exc()))


rule__run_ariba_mlst(
    snakemake.input,
    snakemake.output,
    snakemake.params.sampleComponentObj,
    snakemake.log)
