from bifrostlib import datahandling

def script__run_mlst(reads, sample_file_name, component_file_name, log):
    try:
        db_sample = datahandling.load_sample(sample_file_name)
        db_component = datahandling.load_component(component_file_name)
        log_out = str(log.out_file)
        log_err = str(log.err_file)

        datahandling.log(log_out, "Started {}\n".format("script__run_mlst"))
        species = db_sample["properties"]["species"]

        mlst_species = []
        if species not in db_component["mlst_species_mapping"]:
            datahandling.log(log_out, "cge mlst species: {}\n".format(species))
            shell("touch no_mlst_species_DB")
        else:
            mlst_database_path = db_component["mlst_database_path"]
            mlst_species = db_component["mlst_species_mapping"][species]
            for mlst_entry in mlst_species:
                mlst_entry_path = component + "/" + mlst_entry
                datahandling.log(log_out, "mlst {} on species: {}\n".format(mlst_entry, species))
                shell("if [ -d \"{}\" ]; then rm -r {}; fi".format(mlst_entry_path, mlst_entry_path))
                shell("mkdir {}; mlst.py -x -matrix -s {} -p {} -mp kma -i {} {} -o {} 1> {} 2> {}".format(mlst_entry_path, mlst_entry, mlst_database_path, input.reads[0], input.reads[1], mlst_entry_path, log.out_file, log.err_file))
        shell("touch {}".format(output.complete))
        datahandling.log(log_out, "Done {}\n".format("script__run_mlst"))
    except Exception as e:
        datahandling.log(log_err, str(traceback.format_exc()))
    return 0

script__run_mlst(snakemake.input.reads, snakemake.params.sample_file_name, snakemake.params.component_file_name, snakemake.log)
