#!/usr/bin/env python3
import gzip
import ruamel.yaml
import statistics
import re
import pandas
import os
import vcf
import pkg_resources
import configparser #for checking miseq/nextseq samplesheets
import datetime
import hashlib
config_file = pkg_resources.resource_filename(__name__, "../config/config.yaml")

yaml = ruamel.yaml.YAML(typ='safe')
yaml.default_flow_style = False
with open(config_file, "r") as yaml_stream:
    config = yaml.load(yaml_stream)


def check__sample_sheet():
    #samplesheet = configparser.

    return 0


def md5sum(file):
    with open(file, 'rb') as fh:
        md5sum = hashlib.md5()
        for data in iter(lambda: fh.read(4096), b""):
            md5sum.update(data)
    return md5sum.hexdigest()


def check__run_folder(run_folder, run_info_yaml="run.yaml"):
    run_info = {'samples': {}}
    for file in sorted(os.listdir(run_folder)):
        result = re.search(config["serum"]["read_pattern"], file)
        if result and os.path.isfile(os.path.join(run_folder, file)):
            if result.group("sample_name") not in run_info['samples']:
                run_info['samples'][str(result.group("sample_name"))] = {"count": 0}
            run_info['samples'][result.group("sample_name")]["count"] += 1
            run_info['samples'][result.group("sample_name")][result.group("paired_read_number")] = os.path.join(run_folder, file)
            run_info['samples'][result.group("sample_name")][result.group("paired_read_number") + '_md5sum'] = md5sum(os.path.join(run_folder, file))
    for sample in run_info['samples']:
        if sample not in config["serum"]["samples_to_ignore"]:
            if run_info['samples'][sample]["count"] > 2:
                print("ERROR: at least one sample has more than 2 read_files associted with it, check config['serum']['read_pattern'] for read_pattern")
                return 1
            if run_info['samples'][sample]["count"] == 1:
                print("{} SE read".format(sample))
                return 1
            else:
                # print("{} PE read".format(sample))
                pass
    # return run_info
    with open(run_info_yaml, "w") as output:
        yaml.dump(run_info, output)
        return 0


def check__combine_sample_sheet_with_run_info(sample_sheet_xlsx, run_info_yaml="run.yaml", updated_run_info_yaml="run.yaml"):
    with open(run_info_yaml, "r") as yaml_stream:
        run_info = yaml.load(yaml_stream)
    mapped_columns = {}
    for item in config["serum"]["samplesheet_column_mapping"]:
        mapped_columns[config["serum"]["samplesheet_column_mapping"][item]] = item
    if not os.path.isfile(sample_sheet_xlsx):
        with open(updated_run_info_yaml, "w") as output:
            yaml.dump(run_info, output)
        return 0
    sample_sheet = pandas.read_excel(sample_sheet_xlsx)
    for i, row in sample_sheet.iterrows():
        if config["serum"]["samplesheet_column_mapping"]["sample_name"] not in sample_sheet:
            print("ERROR: sample sheet does not have column")
            return 1
        sample_name = row[str(config["serum"]["samplesheet_column_mapping"]["sample_name"])]
        if sample_name not in run_info['samples']:
            run_info['samples'][sample_name] = {}
        for column in sample_sheet:
            if isinstance(row[column], datetime.date):
                value = row[column].isoformat()
            else:
                value = row[column]
            if column not in mapped_columns:
                run_info['samples'][sample_name][column] = value
            else:
                run_info['samples'][sample_name][mapped_columns[column]] = value
    with open(updated_run_info_yaml, "w") as output:
        yaml.dump(run_info, output)
    return 0


def initialize__run_from_run_info(updated_run_info_yaml="run.yaml", run_status="run_status.csv"):
    with open(updated_run_info_yaml, "r") as yaml_stream:
        run_info = yaml.load(yaml_stream)
    for sample in run_info["samples"]:
        os.makedirs(sample)
        with open(os.path.join(sample, "sample.yaml"), "w") as sample_yaml:
            if sample in config['serum']['samples_to_ignore']:
                run_info["samples"][sample]["status"] = "skipped"
            elif "R1" not in run_info["samples"][sample] or "R2" not in run_info["samples"][sample]:
                run_info["samples"][sample]["status"] = "no reads"
            else:
                run_info["samples"][sample]["status"] = "initialized"
                with open(os.path.join(sample, "cmd"), "w") as command:
                    command.write("snakemake -s ~/code/serumqc/snakefiles/serumqc.snake --config R1_reads={} R2_reads={}".format(run_info["samples"][sample]["R1"], run_info["samples"][sample]["R2"]))
            yaml.dump({"sample": run_info["samples"][sample]}, sample_yaml)
    convert_run_to_status()
    return 0


def initialize_complete():
    with open("init_complete", "w"):
        pass


def convert_run_to_status(run_dir=".", run_status="run_status.csv"):
    run_info = {}
    for directory in os.listdir(run_dir):
        if os.path.isfile(os.path.join(directory, "sample.yaml")):
            with open(os.path.join(directory, "sample.yaml")) as yaml_stream:
                run_info[directory] = yaml.load(yaml_stream)['sample']
    df = pandas.DataFrame.from_dict(run_info, orient='index')
    df.to_csv(run_status)


def start_initialized_samples(run_dir=".", run_status="run_status.csv"):
    for directory in os.listdir(run_dir):
        if os.path.isfile(os.path.join(directory, "sample.yaml")):
            with open(os.path.join(directory, "sample.yaml"), "r") as yaml_stream:
                sample_info = yaml.load(yaml_stream)
            if sample_info['sample']['status'] == 'initialized':
                sample_info['sample']['status'] = 'starting'
                print("running sample")
            with open(os.path.join(directory, "sample.yaml"), "w") as yaml_stream:
                yaml.dump(sample_info, yaml_stream)
    convert_run_to_status()
    return 0


def check__robot_sample_sheet(sample_sheet_xlsx, corrected_sample_sheet_xlsx):
    # change logic to apply to N_WGS and then used N_WGS to replace H_WGS
    df = pandas.read_excel(sample_sheet_xlsx)
    item_rename_dict = {}
    badly_named_samples = df[df["SampleID"].str.contains("^[a-zA-Z0-9\-_]+$") == False]  # samples which fail this have inappropriate characters

    for item in badly_named_samples["SampleID"].tolist():
        print("Renaming {} to {}".format(item, re.sub(r'[^a-zA-Z0-9\-_]', "_", str(item))))

    for item in df["SampleID"].tolist():
        item_rename_dict[item] = re.sub(r'[^a-zA-Z0-9\-_]', "_", str(item))

    df["SampleID"] = df["SampleID"].map(item_rename_dict)

    duplicates = {}
    for item in df["SampleID"].tolist():
        if df["SampleID"].tolist().count(item) > 1:
            print("Duplicate SampleID's exist {}".format(item))
            duplicates[item] = df["SampleID"].tolist().count(item)
            # duplicates have to be done on id basis

    ridiculous_count = 0
    for i, row in df.iterrows():
        if df.loc[i, "SampleID"] in duplicates:
            if df.loc[i, "SampleID"] + "_RENAMED-" + str(duplicates[df.loc[i, "SampleID"]]) in df["SampleID"].tolist():
                ridiculous_count += 1
                df.loc[i, "SampleID"] = "VERY_BAD_SAMPLE_NAME_" + str(ridiculous_count)
            else:
                df.loc[i, "SampleID"] = df.loc[i, "SampleID"] + "_RENAMED-" + str(duplicates[df.loc[i, "SampleID"]])
            duplicates[row["SampleID"]] -= 1
    df.to_excel(corrected_sample_sheet_xlsx)
    return 0


def script__generate_sample_values_from_robot_sample_sheet(sample_sheet_xlsx):

    return 0


def script__kmer_reads(read_files):

    return 0


def script__link_read_files_to_cwd(read_files, check_file):
    os.symlink(read_files[0], os.path.basename(read_files[0]))
    os.symlink(read_files[1], os.path.basename(read_files[1]))

    with open(check_file, "w") as output:
        output.write("{} -> {}\n".format(read_files[0], os.path.basename(read_files[0])))
        output.write("{} -> {}\n".format(read_files[1], os.path.basename(read_files[1])))

    return 0


def check__minimum_reads(read_files, minimum_read_check_file):
    """
    read_files = Raw read files
    """
    minimum_number_of_reads = config["serum"]["minimum_number_of_reads"]
    counts = [0] * len(read_files)
    for i, file in enumerate(read_files):
        with gzip.open(file) as read_file:
            for j, line in enumerate(read_file):
                pass
            counts[i] = j + 1

    if (counts[0] == counts[1] and
        counts[0] % 4 == 0 and
        counts[0] / 4 > minimum_number_of_reads / 2):
        with open(minimum_read_check_file, "w") as output:
            output.write("{} total reads and matching number of reads\n".format(counts[0] / 4 * 2))

    return 0


def script__summarize_read_stats(read_files, summarize_read_stats_yaml):
    """
    Needs testing

    read_files is assumed to come as a list of tuples containing:
    input reads (R1,R2),
    trimmed reads (R1,R2),
    normalized reads (R1,R2)
    """
    read_info_dict = {"raw": {}, "trimmed": {}, "normalized": {}}
    read_info_index = ["raw", "trimmed", "normalized"]

    for i, read_pair in enumerate(read_files):
        for j, read_file in enumerate(read_pair):
            name = "R1"
            if j % 2 == 0:
                name = "R2"
            read_lengths = []
            with gzip.open(read_file) as reads:
                for k, line in enumerate(reads):
                    if k % 4 == 3:
                        read_lengths.append(len(line))
            read_file_dict = {}
            read_file_dict[os.path.basename(name)] = {
                "file": os.path.basename(read_file),
                "reads": len(read_lengths),
                "mean": statistics.mean(read_lengths),
                "mode": statistics.mode(read_lengths),
                "min": min(read_lengths),
                "max": max(read_lengths)
            }
            read_info_dict[read_info_index[i]].update(read_file_dict)

    output_dict = {"read_stats": read_info_dict}
    with open(summarize_read_stats_yaml, "w") as output_file:
        yaml.dump(output_dict, output_file)

    return 0


def script__summarize_kraken_report(kraken_report_file, summarize_contamination_reads_yaml):
    """
    summarize kraken report
    """
    cutoff_threshold = config["kraken"]["summary"]["cutoff_threshold"]
    percent_unclassified = 0.0
    detected_species = "NA"

    kraken_contents = []
    with open(kraken_report_file, "r") as kraken_report:
        for line in kraken_report:
            kraken_contents.append(line.strip().split("\t"))

        # formatting
    previous_spaces = 0
    for i, item in enumerate(kraken_contents):
        item[0] = float(item[0].strip())
        item.append(len(item[5]) - len(item[5].lstrip(' ')))
        ########
        if item[0] < cutoff_threshold:
            item.append("below threshold")
        elif item[6] > previous_spaces and previous_spaces != 0:
            item.append("below level")
        elif item[6] < previous_spaces:
            item.append("OK")
            previous_spaces = 0
        elif item[3] == "S":
            item.append("OK")
            previous_spaces = item[6]
        else:
            item.append("OK")

    df = pandas.DataFrame(kraken_contents, columns=["classified", "total_reads", "reads_at_level", "level", "ncbi", "name", "spaces", "include_in_summary"])

    # Unclassified is always the first position
    percent_unclassified = float(df[df["level"] == "U"].loc[0]["classified"])

    filtered_df = df[df["include_in_summary"] == "OK"]

    # Retrieve dominant species if present
    detected_species_index = -1
    if not filtered_df[filtered_df["level"] == "S"].empty:
        detected_species_index = filtered_df[filtered_df["level"] == "S"].index[0]
        detected_species = filtered_df[filtered_df["level"] == "S"].loc[detected_species_index]["name"].lstrip(' ')

    # Retrieve distinct taxa groups over threshold
    taxa_detected = []
    previous_spaces = filtered_df.loc[detected_species_index - 1]["spaces"]
    previous_index = detected_species_index - 1
    if detected_species_index >= 0:
        for index, row in filtered_df.loc[detected_species_index - 1:].iterrows():
            if row["level"] == "S":
                pass
            elif row["spaces"] < previous_spaces:
                taxa_detected.append(filtered_df.loc[previous_index]["name"].lstrip(' ') + " " + filtered_df.loc[previous_index]["level"])
                previous_spaces = row["spaces"]
                previous_index = index
            else:
                previous_spaces = row["spaces"]
                previous_index = index
        # One item left that needs to be added, will always happen if a species is detected
        taxa_detected.append(filtered_df.loc[previous_index]["name"].lstrip(' ') + " " + filtered_df.loc[previous_index]["level"])

    output_dict = {
        "contamination_reads": {
            "percent_unclassified": percent_unclassified,
            "species_detected": detected_species,
            "taxonomies_detected": taxa_detected,
            "number_of_taxa_in_sample": len(taxa_detected)
        }
    }
    with open(summarize_contamination_reads_yaml, "w") as output_file:
        yaml.dump(output_dict, output_file)

    return 0


def script__summarize_depth(depth_file, summarize_contig_depth_yaml, summarize_binned_depth_yaml):
    """

    """
    depth_dict = {}
    with open(depth_file, "r") as input_file:
        for line in input_file:
            contig = line.split("\t")[0]
            depth = int(line.split("\t")[2].strip())
            if contig not in depth_dict:
                depth_dict[contig] = {}
            if depth in depth_dict[contig]:
                depth_dict[contig][depth] += 1
            else:
                depth_dict[contig][depth] = 1

    contig_depth_summary_dict = {}
    for contig in depth_dict:
        total_depth = 0
        total_length = 0
        for depth in depth_dict[contig]:
            length = depth_dict[contig][depth]
            total_depth = total_depth + (length * depth)
            total_length = total_length + length
        contig_depth_summary_dict[contig] = {
            "coverage": total_depth / total_length,
            "total_depth": total_depth,
            "total_length": total_length
        }

    # Removing as it's also looped in vcf so no need to check twice
    # dict is made now to cycle over on range and on contigs
    binned_depth_summary_dict = {}
    depth_limits = config["serum"]["summarize"]["depth_range"]
    depth_range = list(range(depth_limits[0], depth_limits[1]))
    for bound in depth_range:
        binned_depth_summary_dict[bound] = 0

    for contig in depth_dict:
        for depth in depth_dict[contig]:
                limit = depth
                if depth > depth_limits[1]:
                    limit = depth_limits[1]
                for i in range(depth_limits[0], limit):
                    if depth >= i:
                        binned_depth_summary_dict[i] += depth_dict[contig][depth]

    output_dict = {"contig_depth": contig_depth_summary_dict}
    with open(summarize_contig_depth_yaml, "w") as output_file:
        yaml.dump(output_dict, output_file)

    output_dict = {"binned_depth": binned_depth_summary_dict}
    with open(summarize_binned_depth_yaml, "w") as output_file:
        yaml.dump(output_dict, output_file)

    return 0


def check__detect_species(summarize_kraken_report_yaml, species_file):
    """
    Refactor this away
    """
    with open(summarize_kraken_report_yaml, "r") as yaml_stream:
        serum_db = yaml.load(yaml_stream)
    with open(species_file, "w") as outfile:
        outfile.write(serum_db["contamination_reads"]["species_detected"] + "\n")

    return 0


def script__summarize_spades_log(spades_log_file, summarize_insert_deviation_yaml):
    insert_deviation_dict = {}
    with open(spades_log_file) as spades_log:
        content = spades_log.read()
        insert_deviation_dict["insert_size"] = re.search("Insert size = (.*?), deviation = (.*?),", content, re.MULTILINE).group(1)
        insert_deviation_dict["deviation"] = re.search("Insert size = (.*?), deviation = (.*?),", content, re.MULTILINE).group(2)

    output_dict = {"insert_deviation": insert_deviation_dict}
    with open(summarize_insert_deviation_yaml, "w") as output_file:
        yaml.dump(output_dict, output_file)
    return 0


def script__summarize_mlst(mlst_file, serum_db):
    output_dict = {"mlst": ""}
    with open(serum_db, "w") as output_file:
        yaml.dump(output_dict, output_file)

    return 0


def script__summarize_quast_report(quast_report_file, summarize_assembly_yaml):
    assembly_stats_dict = {}
    with open(quast_report_file, "r") as quast_report:
        content = quast_report.read()
        assembly_stats_dict["N50"] = re.search("N50\t(?P<N50>[0-9]+)", content).group(1)
        assembly_stats_dict["N75"] = re.search("N75\t(?P<N75>[0-9]+)", content).group(1)
        assembly_stats_dict["L50"] = re.search("L50\t(?P<L50>[0-9]+)", content).group(1)
        assembly_stats_dict["L75"] = re.search("L75\t(?P<L75>[0-9]+)", content).group(1)

    output_dict = {"assembly": assembly_stats_dict}
    with open(summarize_assembly_yaml, "w") as output_file:
        yaml.dump(output_dict, output_file)

    return 0


def script__summarize_pilon_vcf(pilon_vcf_file, summarize_ambiguous_snp_yaml):
    ambiguous_snp_dict = {}

    pilon_vcf = vcf.Reader(open(pilon_vcf_file, "r"))
    depth_limits = config["serum"]["summarize"]["depth_range"]
    depth_range = list(range(depth_limits[0], depth_limits[1]))
    frequency_range = [round(i / 100, 2) for i in range(0, 51)]
    for bound in depth_range:
        ambiguous_snp_dict[bound] = {}
        for frequency in frequency_range:
            ambiguous_snp_dict[bound][frequency] = 0

    vcf_dict = {}
    for i, record in enumerate(pilon_vcf):
        record_dict = {}
        if not ("DP" in record.INFO and "AF" in record.INFO):
            pass
        else:
            record_dict["depth"] = int(record.INFO["DP"])
            if round(sum(record.INFO["AF"]), 2) < 0.50:
                record_dict["AF"] = round(sum(record.INFO["AF"]), 2)
            else:
                record_dict["AF"] = round(1.00 - sum(record.INFO["AF"]), 2)

            if record.affected_end - record.affected_start != 1:
                record_dict["indel_or_deletion"] = True
            else:
                record_dict["indel_or_deletion"] = False
            vcf_dict[i] = record_dict

    for vcy_key, vcf_entry in vcf_dict.items():
        if not vcf_entry["indel_or_deletion"]:
            depth_max = vcf_entry["depth"]
            if depth_max > depth_limits[1]:
                depth_max = depth_limits[1]
            for depth in range(depth_limits[0], depth_max):
                for allele_frequency in frequency_range:
                    if vcf_entry["AF"] > allele_frequency:
                        ambiguous_snp_dict[depth][allele_frequency] += 1

    output_dict = {"ambiguous_snp": ambiguous_snp_dict}
    with open(summarize_ambiguous_snp_yaml, "w") as output_file:
        yaml.dump(output_dict, output_file)

    return 0


def script__summarize_kraken_contigs(kraken_contigs_file, summarize_contamination_contigs_yaml):
    kraken_name_file = config["kraken"]["names"]
    name_dict = {'0': "unclassified"}
    with open(kraken_name_file, "r") as name_file:
        for line in name_file:
            info = line.split("|")
            for i in range(0, len(info)):
                info[i] = info[i].strip()
            if info[0] not in name_dict:
                name_dict[info[0]] = info[1]
            if info[3] == "scientific name":
                name_dict[info[0]] = info[1]

    contaminations_contigs_dict = {}
    with open(kraken_contigs_file, "r") as kraken_file:
        for line in kraken_file:
            info = line.split("\t")
            contaminations_contigs_dict[info[1]] = name_dict[info[2]]

    output_dict = {"contamination_contigs": contaminations_contigs_dict}
    with open(summarize_contamination_contigs_yaml, "w") as output_file:
        yaml.dump(output_dict, output_file)

    return 0


def script__summarize_summaries_to_serum(yaml_files, serumqc_summary_yaml):
    """
    """
    serumqc_summary = {}
    for file in yaml_files:
        with open(file, "r") as yaml_stream:
            partial_info = (yaml.load(yaml_stream))
            for key in partial_info:
                if key in serumqc_summary:
                    print("Yaml key error, serumqc key: {}".format(key))
                serumqc_summary[key] = partial_info[key]

    with open(serumqc_summary_yaml, "w") as output_file:
        yaml.dump(serumqc_summary, output_file)
    # combine yaml files
    return 0


def qc_yaml(serumqc_summary_yaml, serumqc_yaml):
    """

    """
    serumqc_dict = {}

    with open(serumqc_summary_yaml, "r") as yaml_stream:
        serumqc_summary = yaml.load(yaml_stream)
    serumqc_summary["read_stats"]["raw"]["R1"]
    serumqc_summary["read_stats"]["raw"]["R2"]

    # Number of reads
    serumqc_dict["num_of_reads"] = serumqc_summary["read_stats"]["raw"]["R1"]["reads"] + serumqc_summary["read_stats"]["raw"]["R2"]["reads"]
    serumqc_dict["read_length"] = serumqc_summary["read_stats"]["raw"]["R1"]["mode"]
    serumqc_dict["species_detected"] = serumqc_summary["contamination_reads"]["species_detected"]

    serumqc_dict["contigs"] += 0
    total_depth = 0
    total_length = 0
    for contig in serumqc_summary["contig_depth"]:
        if contig["total_length"] >= config["qc"]["min_length"]:
            serumqc_dict["contigs"] += 1
            total_depth += contig["total_depth"]
            total_length += contig["total_length"]

    serumqc_dict["length"] = total_length
    serumqc_dict["coverage"] = total_depth / total_length
    serumqc_dict["N50"] = serumqc_summary["assembly"]["N50"]
    serumqc_dict["N75"] = serumqc_summary["assembly"]["N75"]

    base_depth = serumqc_summary["binned_depth"][config["serum"]["qc"]["base_depth"]]
    min_depth = serumqc_summary["binned_depth"][config["serum"]["qc"]["min_depth"]]
    serumqc_summary["binned_depth"][config["serum"]["qc"]["recommended_depth"]]

    # sample_name
    # status?
    # supplying lab
    # emails (initials)
    # run_name
    # qc action
    # output directory
    # number of reads
    # read length
    # detected species ncbi
    # mlst type
    # mlst alleles
    # coverage base
    # coverage compare
    # length
    # contigs
    # coverage
    # N50
    # N75
    # plate name
    # comments

    return 0

def format_contigs():
    # todo
    return 0