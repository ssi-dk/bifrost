import os
from bifrostlib import datahandling


def extract_contigs_sum_cov(sampleComponentObj):
    summary, results = sampleComponentObj.get_summary_and_results()
    options = sampleComponentObj.get_options()
    file_path = os.path.join(sampleComponentObj.get_component_name(), "contigs.sum.cov")
    key = file_path.replace(".","_").replace("$","_")
    results[key] = {}

    yaml = datahandling.load_yaml(file_path)
    for bin_value in options["cov_bin_values"]:
        total_length = 0
        total_depth = 0
        total_contigs = 0
        for contig in yaml["contig_depth"]:
            if yaml["contig_depth"][contig]["coverage"] >= float(bin_value):
                total_length += yaml["contig_depth"][contig]["total_length"]
                total_depth += yaml["contig_depth"][contig]["total_depth"]
                total_contigs += 1
        results[key]["bin_contigs_at_{}x".format(bin_value)] = total_contigs
        results[key]["bin_length_at_{}x".format(bin_value)] = total_length
        results[key]["bin_coverage_at_{}x".format(bin_value)] = float(total_depth / total_length)
        summary["bin_contigs_at_{}x".format(bin_value)] = total_contigs
        summary["bin_length_at_{}x".format(bin_value)] = total_length
        summary["bin_coverage_at_{}x".format(bin_value)] = float(total_depth / total_length)
    return (summary, results)


def extract_contigs_bin_cov(db, file_path, key, temp_data):
    summary, results = sampleComponentObj.get_summary_and_results()
    options = sampleComponentObj.get_options()
    file_path = os.path.join(sampleComponentObj.get_component_name(), "contigs.bin.cov")
    key = file_path.replace(".", "_").replace("$", "_")
    results[key] = {}

    yaml = datahandling.load_yaml(file_path)
    results[key] = yaml
    for bin_value in GLOBAL_cov_bin_values:
        summary["raw_length_at_{}x".format(bin_value)] = yaml["binned_depth"][bin_value - 1]
    return (summary, results)


def extract_bbuk_log(db, file_path, key, temp_data):
    import json
    buffer = datahandling.read_buffer(file_path)
    bbduk_dict = json.loads("{" + buffer.partition("{")[2])
    db["results"][key] = bbduk_dict
    db["results"][key]["filtered_reads_num"] = db["results"][key]["readsIn"] - db["results"][key]["readsRemoved"]
    db["summary"]["filtered_reads_num"] = db["results"][key]["filtered_reads_num"]
    return db


def extract_quast_report(db, file_path, key, temp_data):
    buffer = datahandling.read_buffer(file_path)
    db["results"][key]["GC"] = float(re.search("GC \(%\)\t([0-9]+[\.]?[0-9]*)", buffer, re.MULTILINE).group(1))
    db["results"][key]["N50"] = int(re.search("N50\t([0-9]+)", buffer, re.MULTILINE).group(1))
    db["results"][key]["N75"] = int(re.search("N75\t([0-9]+)", buffer, re.MULTILINE).group(1))
    db["results"][key]["L50"] = int(re.search("L50\t([0-9]+)", buffer, re.MULTILINE).group(1))
    db["results"][key]["L75"] = int(re.search("L75\t([0-9]+)", buffer, re.MULTILINE).group(1))
    db["summary"]["GC"] = db["results"][key]["GC"]
    db["summary"]["N50"] = db["results"][key]["N50"]
    return db


def extract_contig_variants(db, file_path, key, temp_data):
    yaml = datahandling.load_yaml(file_path)
    db["results"][key] = yaml
    db["summary"]["snp_filter_10x_10%"] = yaml["variant_table"][9][9]
    db["summary"]["snp_filter_indels"] = yaml["indels"]
    db["summary"]["snp_filter_deletions"] = yaml["deletions"]
    return db


def extract_contig_stats(db, file_path, key, temp_data):
    buffer = datahandling.read_buffer(file_path)
    for line in buffer.split("\n"):
        if line.startswith("SN"):
            db["results"][key][re.sub('[^A-Za-z0-9\s]+', '', line.split("\t")[1]).replace(" ", "_").rstrip("_")] = line.split("\t")[2]
    db["summary"]["raw_total_sequences"] = db["results"][key]["raw_total_sequences"]
    db["summary"]["reads_mapped"] = db["results"][key]["reads_mapped"]
    db["summary"]["reads_unmapped"] = db["results"][key]["reads_unmapped"]
    db["summary"]["insert_size_average"] = db["results"][key]["insert_size_average"]
    db["summary"]["insert_size_standard_deviation"] = db["results"][key]["insert_size_standard_deviation"]
    return db


def script__datadump(output, sample_file, component_file, sample_component_file, log):
    try:
        output = str(output)
        log_out = str(log.out_file)
        log_err = str(log.err_file)
        db_sample = datahandling.load_sample(sample_file)
        db_component = datahandling.load_component(component_file)
        db_sample_component = datahandling.load_sample_component(sample_component_file)
        this_function_name = sys._getframe().f_code.co_name
        global GLOBAL_component_name
        GLOBAL_component_name = db_component["name"]
        global GLOBAL_cov_bin_values
        GLOBAL_cov_bin_values = db_component["options"]["cov_bin_values"]
        datahandling.write_log(log_out, "Started {}\n".format(this_function_name))

        # Save files to DB
        datahandling.save_files_to_db(db_component["db_values_changes"]["files"], sample_component_id=db_sample_component["_id"])

        # Initialization of values, summary and reporter are also saved into the sample
        db_sample_component["summary"] = {"component": {"_id": db_component["_id"], "_date": datetime.datetime.utcnow()}}
        db_sample_component["results"] = {}
        db_sample_component["reporter"] = {}

        db_sample_component = datahandling.datadump_template(extract_contigs_sum_cov, db_sample_component, file_path=os.path.join(GLOBAL_component_name, "contigs.sum.cov"))
        db_sample_component = datahandling.datadump_template(extract_contigs_bin_cov, db_sample_component, file_path=os.path.join(GLOBAL_component_name, "contigs.bin.cov"))
        db_sample_component = datahandling.datadump_template(extract_bbuk_log, db_sample_component, file_path=os.path.join(GLOBAL_component_name, "log/setup__filter_reads_with_bbduk.err.log"), )
        db_sample_component = datahandling.datadump_template(extract_quast_report, db_sample_component, file_path=os.path.join(GLOBAL_component_name, "quast/report.tsv"))
        db_sample_component = datahandling.datadump_template(extract_contig_variants, db_sample_component, file_path=os.path.join(GLOBAL_component_name, "contigs.variants"))
        db_sample_component = datahandling.datadump_template(extract_contig_stats, db_sample_component, file_path=os.path.join(GLOBAL_component_name, "contigs.stats"))

        datahandling.save_sample_component_to_file(db_sample_component, sample_component_file)
        db_sample["properties"]["denovo_assembly"] = db_sample_component["summary"]
        datahandling.save_sample_to_file(db_sample, sample_file)
        open(output, 'w+').close()  # touch file

    except Exception:
        datahandling.write_log(log_out, "Exception in {}\n".format(this_function_name))
        datahandling.write_log(log_err, str(traceback.format_exc()))
        raise Exception
        return 1

    finally:
        datahandling.write_log(log_out, "Done {}\n".format(this_function_name))
        return 0


def datadump(sampleComponentObj, log):
    sampleComponentObj.start_data_dump(log=log)
    sampleComponentObj.run_data_dump_on_function(extract_contigs_sum_cov, log=log)
    sampleComponentObj.run_data_dump_on_function(extract_contigs_bin_cov, log=log)
    sampleComponentObj.run_data_dump_on_function(extract_bbuk_log, log=log)
    sampleComponentObj.run_data_dump_on_function(extract_quast_report, log=log)
    sampleComponentObj.run_data_dump_on_function(extract_contig_variants, log=log)
    sampleComponentObj.run_data_dump_on_function(extract_contig_stats, log=log)
    sampleComponentObj.end_data_dump(log=log)


datadump(
    snakemake.params.sampleComponentObj,
    snakemake.log)
