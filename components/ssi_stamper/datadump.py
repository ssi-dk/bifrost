import pkg_resources
import datetime
import os
import re
import sys
import traceback
from bifrostlib import datahandling




def test__sample__has_reads_files(db, file_path, key, temp_data):
    try:
        sample_db = temp_data["sample_db"]
        component_db = temp_data["component_db"]
        this_function_name = sys._getframe().f_code.co_name

        test = {
            "name": this_function_name,
            "display_name": "No reads",
            "effect": "core facility",
            "value": "",
            "status": "",
            "reason": ""
        }
        if sample_db["reads"]["R1"] == "":
            test["status"] = "fail"
            test["reason"] = "Read path is empty"
        else:
            test["status"] = "pass"
            test["reason"] = ""

    except KeyError as e:
        test["status"] = "fail"
        test["reason"] = "Database KeyError {} in function {}: ".format(e.args[0], this_function_name)

    finally:
        db["results"][this_function_name] = test
        return db


def test__species_detection__main_species_level(db, file_path, key, temp_data):
    try:
        sample_db = temp_data["sample_db"]
        component_db = temp_data["component_db"]
        this_function_name = sys._getframe().f_code.co_name

        test = {
            "name": this_function_name,
            "display_name": "Multiple species detected",
            "effect": "supplying lab",
            "value": "",
            "status": "",
            "reason": ""
        }
        test["value"] = round(sample_db["properties"]["species_detection"]["percent_classified_species_1"] + sample_db["properties"]["species_detection"]["percent_unclassified"], 3)
        if test["value"] < component_db["options"]["min_species"]:
            test["status"] = "fail"
            test["reason"] = "Value ({}) is below threshold ({})".format(test["value"], component_db["options"]["min_species"])
        else:
            test["status"] = "pass"
            test["reason"] = ""

    except KeyError as e:
        test["status"] = "fail"
        test["reason"] = "Database KeyError {} in function {}: ".format(e.args[0], this_function_name)

    finally:
        db["results"][this_function_name] = test
        return db


def test__species_detection__unclassified_level(db, file_path, key, temp_data):
    try:
        sample_db = temp_data["sample_db"]
        component_db = temp_data["component_db"]
        this_function_name = sys._getframe().f_code.co_name

        test = {
            "name": this_function_name,
            "display_name": "High unclassified",
            "effect": "supplying lab",
            "value": "",
            "status": "",
            "reason": ""
        }
        test["value"] = round(sample_db["properties"]["species_detection"]["percent_unclassified"], 3)
        if test["value"] >= component_db["options"]["max_unclassified"]:
            test["status"] = "fail"
            test["reason"] = "Value ({}) is above threshold ({})".format(test["value"], component_db["options"]["max_unclassified"])
        else:
            test["status"] = "pass"
            test["reason"] = ""

    except KeyError as e:
        test["status"] = "fail"
        test["reason"] = "Database KeyError {} in function {}: ".format(e.args[0], this_function_name)

    finally:
        db["results"][this_function_name] = test
        return db


def test__component__species_in_db(db, file_path, key, temp_data):
    try:
        sample_db = temp_data["sample_db"]
        component_db = temp_data["component_db"]
        this_function_name = sys._getframe().f_code.co_name

        test = {
            "name": this_function_name,
            "display_name": "No species submitted - using default values",
            "effect": "supplying lab",
            "value": "",
            "status": "",
            "reason": ""
        }
        test["value"] = sample_db["properties"].get("species", None)

        if test["value"] not in component_db["options"]["species_qc_value_mapping"]:
            test["status"] = "fail"
            test["reason"] = "Detected species not in bifrost db. Can't estimate proper QC values."
        else:
            test["status"] = "pass"
            test["reason"] = ""

    except KeyError as e:
        test["status"] = "fail"
        test["reason"] = "Database KeyError {} in function {}: ".format(e.args[0], this_function_name)

    finally:
        db["results"][this_function_name] = test
        return db


def test__sample__species_provided_is_detected(db, file_path, key, temp_data):
    try:
        sample_db = temp_data["sample_db"]
        component_db = temp_data["component_db"]
        this_function_name = sys._getframe().f_code.co_name

        test = {
            "name": this_function_name,
            "display_name": "Detected species mismatch",
            "effect": "supplying lab",
            "value": "",
            "status": "",
            "reason": ""
        }
        test["value"] = sample_db["properties"].get("provided_species", None)
        species = temp_data["sample_db"]["properties"]["species"]
        if test["value"] is None:
            test["status"] = "pass"
            test["reason"] = "No submitted species"
        elif test["value"] not in component_db["options"]["species_qc_value_mapping"]:
            test["status"] = "pass"
            test["reason"] = "Submitted species not in db"
        elif species != test["value"]:
            test["status"] = "fail"
            test["reason"] = "Detected species ({}) different than expected ({})".format(species, test["value"])
        else:
            test["status"] = "pass"
            test["reason"] = ""

    except KeyError as e:
        test["status"] = "fail"
        test["reason"] = "Database KeyError {} in function {}: ".format(e.args[0], this_function_name)

    finally:
        db["results"][this_function_name] = test
        return db


def test__denovo_assembly__genome_size_at_1x(db, file_path, key, temp_data):
    try:
        sample_db = temp_data["sample_db"]
        component_db = temp_data["component_db"]
        this_function_name = sys._getframe().f_code.co_name

        test = {
            "name": this_function_name,
            "display_name": "Atypical genome size (1x)",
            "effect": "supplying lab",
            "value": "",
            "status": "",
            "reason": ""
        }
        test["value"] = sample_db["properties"]["denovo_assembly"]["bin_contigs_at_1x"]
        species = temp_data["sample_db"]["properties"]["species"]
        if species not in component_db["options"]["species_qc_value_mapping"]:
            species = "default"
        min_length = component_db["options"]["species_qc_value_mapping"][species]["min_length"]
        max_length = component_db["options"]["species_qc_value_mapping"][species]["max_length"]
        if min_length < test["value"] < max_length:
            test["status"] = "pass"
            test["reason"] = ""
        else:
            test["status"] = "fail"
            test["reason"] = "Value ({}) below or above expected ({}, {})".format(test["value"], min_length, max_length)

    except KeyError as e:
        test["status"] = "fail"
        test["reason"] = "Database KeyError {} in function {}: ".format(e.args[0], this_function_name)

    finally:
        db["results"][this_function_name] = test
        return db


def test__denovo_assembly__genome_size_at_10x(db, file_path, key, temp_data):
    try:
        sample_db = temp_data["sample_db"]
        component_db = temp_data["component_db"]
        this_function_name = sys._getframe().f_code.co_name

        test = {
            "name": this_function_name,
            "display_name": "Atypical genome size (10x)",
            "effect": "supplying lab",
            "value": "",
            "status": "",
            "reason": ""
        }
        test["value"] = sample_db["properties"]["denovo_assembly"]["bin_contigs_at_10x"]
        species = temp_data["sample_db"]["properties"]["species"]
        if species not in component_db["options"]["species_qc_value_mapping"]:
            species = "default"
        min_length = component_db["options"]["species_qc_value_mapping"][species]["min_length"]
        max_length = component_db["options"]["species_qc_value_mapping"][species]["max_length"]
        if min_length < test["value"] < max_length:
            test["status"] = "pass"
            test["reason"] = ""
        else:
            test["status"] = "fail"
            test["reason"] = "Value ({}) below or above expected ({}, {})".format(test["value"], min_length, max_length)

    except KeyError as e:
        test["status"] = "fail"
        test["reason"] = "Database KeyError {} in function {}: ".format(e.args[0], this_function_name)

    finally:
        db["results"][this_function_name] = test
        return db


def test__denovo_assembly__genome_size_difference_1x_10x(db, file_path, key, temp_data):
    try:
        sample_db = temp_data["sample_db"]
        component_db = temp_data["component_db"]
        this_function_name = sys._getframe().f_code.co_name

        test = {
            "name": this_function_name,
            "display_name": "Atypical genome size difference (1x - 10x)",
            "effect": "supplying lab",
            "value": "",
            "status": "",
            "reason": ""
        }
        test["value"] = sample_db["properties"]["denovo_assembly"]["bin_contigs_at_1x"] - sample_db["properties"]["denovo_assembly"]["bin_contigs_at_10x"]
        max_size_difference = component_db["options"]["max_size_difference_for_1x_and_10x"]
        if test["value"] < max_size_difference:
            test["status"] = "pass"
            test["reason"] = ""
        else:
            test["status"] = "fail"
            test["reason"] = "Value ({}) above expected ({})".format(test["value"], max_size_difference)

    except KeyError as e:
        test["status"] = "fail"
        test["reason"] = "Database KeyError {} in function {}: ".format(e.args[0], this_function_name)

    finally:
        db["results"][this_function_name] = test
        return db


def test__denovo_assembly__genome_average_coverage(db, file_path, key, temp_data):
    try:
        sample_db = temp_data["sample_db"]
        component_db = temp_data["component_db"]
        this_function_name = sys._getframe().f_code.co_name

        test = {
            "name": this_function_name,
            "display_name": "Atypical genome size difference (1x - 10x)",
            "effect": "supplying lab",
            "value": "",
            "status": "",
            "reason": ""
        }
        test["value"] = round(sample_db["properties"]["denovo_assembly"]["bin_coverage_at_1x"], 3)
        average_coverage_fail = component_db["options"]["average_coverage_fail"]
        average_coverage_low = component_db["options"]["average_coverage_low"]
        average_coverage_warn = component_db["options"]["average_coverage_warn"]

        if test["value"] < average_coverage_fail:
            test["status"] = "fail"
            test["reason"] = "Lack of reads ({} < {})".format(test["value"], average_coverage_fail)
        elif test["value"] < average_coverage_low:
            test["status"] = "fail"
            test["reason"] = "Not enough reads ({} < {})".format(test["value"], average_coverage_low)
        elif test["value"] < average_coverage_warn:
            test["status"] = "fail"
            test["reason"] = "Low reads ({} < {})".format(test["value"], average_coverage_warn)
            test["effect"] = "supplying lab"
        else:
            test["status"] = "pass"
            test["reason"] = ""

    except KeyError as e:
        test["status"] = "fail"
        test["reason"] = "Database KeyError {} in function {}: ".format(e.args[0], this_function_name)

    finally:
        db["results"][this_function_name] = test
        return db


def test__denovo_assembly__minimum_read_number(db, file_path, key, temp_data):
    try:
        sample_db = temp_data["sample_db"]
        component_db = temp_data["component_db"]
        this_function_name = sys._getframe().f_code.co_name

        test = {
            "name": this_function_name,
            "display_name": "Number of filtered reads below minimum",
            "effect": "core facility",
            "value": "",
            "status": "",
            "reason": ""
        }
        test["value"] = sample_db["properties"]["denovo_assembly"]["filtered_reads_num"]
        number_of_reads_fail = component_db["options"]["number_of_reads_fail"]

        if test["value"] < number_of_reads_fail:
            test["status"] = "fail"
            test["reason"] = "Filtered reads below minimum ({} < {})".format(num_reads, number_of_reads_fail)
        else:
            test["status"] = "pass"
            test["reason"] = ""

    except KeyError as e:
        test["status"] = "fail"
        test["reason"] = "Database KeyError {} in function {}: ".format(e.args[0], this_function_name)

    finally:
        db["results"][this_function_name] = test
        return db


def evaluate_tests_and_stamp(db, file_path, key, temp_data):
    sample_db = temp_data["sample_db"]
    core_facility = False
    supplying_lab = False
    for test in db["results"]:
        if db["results"][test]["status"] == "fail" or db["results"][test]["status"] == "undefined":
            if db["results"][test]["effect"] == "supplying lab":
                supplying_lab = True
            elif db["results"][test]["effect"] == "core facility":
                core_facility = True
    if (sample_db["properties"]["provided_species"] == sample_db["properties"]["detected_species"] and \
        db["results"]["test__denovo_assembly__genome_average_coverage"]["status"] == "fail" and \
        db["results"]["test__denovo_assembly__genome_average_coverage"]["effect"] == "supplyinh lab" and \
        db["results"]["test__denovo_assembly__genome_size_difference_1x_10x"]["status"] == "fail" and \
        db["results"]["test__sample__species_provided_is_detected"]["status"] == "pass" and \
        db["results"]["test__denovo_assembly__genome_size_at_1x"]["status"] == "pass"):
            core_facility = True
    action = "pass:OK"
    if supplying_lab:
        action = "fail:supplying lab"
    if core_facility:
        action = "fail:core facility"

    db["stamp"] = {
        "name": "ssi_stamper",
        "value": action,
        "date": datetime.datetime.utcnow()
    }
    return db


def generate_summary(db, file_path, key, temp_data):
    for test in db["results"]:
        db["summary"][db["results"][test]["name"]] = "{}:{}:{}".format(db["results"][test]["status"], db["results"][test]["reason"], db["results"][test]["value"])
    return db


def script__datadump(output, sample_file, component_file, sample_component_file, log):
    try:
        output = str(output)
        log_out = str(log.out_file)
        log_err = str(log.err_file)
        sample_db = datahandling.load_sample(sample_file)
        component_db = datahandling.load_component(component_file)
        db_sample_component = datahandling.load_sample_component(sample_component_file)
        this_function_name = sys._getframe().f_code.co_name
        global GLOBAL_component_name
        GLOBAL_component_name = component_db["name"]

        datahandling.write_log(log_out, "Started {}\n".format(this_function_name))

        # Save files to DB
        datahandling.save_files_to_db(component_db["db_values_changes"]["files"], sample_component_id=db_sample_component["_id"])

        # Initialization of values, summary and reporter are also saved into the sample
        db_sample_component["summary"] = {"component": {"_id": component_db["_id"], "_date": datetime.datetime.utcnow()}}
        db_sample_component["results"] = {}
        db_sample_component["reporter"] = {}  # Currently unused, set to dict of component config path when used
#---Unique to component: start----------------------------------------------------------------------
        db_sample_component["tests"] = {}
        sample_db["stamps"] = sample_db.get("stamps", {})
        sample_db["stamps"]["stamp_list"] = sample_db["stamps"].get("stamp_list", [])

        # Variables being used
        working_temp_data = {"sample_db": sample_db, "component_db": component_db}

        db_sample_component = datahandling.datadump_template(test__sample__has_reads_files, db_sample_component, temp_data=working_temp_data)
        db_sample_component = datahandling.datadump_template(test__species_detection__main_species_level, db_sample_component, temp_data=working_temp_data)
        db_sample_component = datahandling.datadump_template(test__species_detection__unclassified_level, db_sample_component, temp_data=working_temp_data)
        db_sample_component = datahandling.datadump_template(test__component__species_in_db, db_sample_component, temp_data=working_temp_data)
        db_sample_component = datahandling.datadump_template(test__sample__species_provided_is_detected, db_sample_component, temp_data=working_temp_data)
        db_sample_component = datahandling.datadump_template(test__denovo_assembly__genome_size_at_1x, db_sample_component, temp_data=working_temp_data)
        db_sample_component = datahandling.datadump_template(test__denovo_assembly__genome_size_at_10x, db_sample_component, temp_data=working_temp_data)
        db_sample_component = datahandling.datadump_template(test__denovo_assembly__genome_size_difference_1x_10x, db_sample_component, temp_data=working_temp_data)
        db_sample_component = datahandling.datadump_template(test__denovo_assembly__genome_average_coverage, db_sample_component, temp_data=working_temp_data)
        db_sample_component = datahandling.datadump_template(test__denovo_assembly__minimum_read_number, db_sample_component, temp_data=working_temp_data)
        db_sample_component = datahandling.datadump_template(evaluate_tests_and_stamp, db_sample_component, temp_data=working_temp_data)
        db_sample_component = datahandling.datadump_template(generate_summary, db_sample_component, temp_data=working_temp_data)

        sample_db["stamps"]["ssi_stamper"] = db_sample_component["stamp"]
        sample_db["stamps"]["stamp_list"].append(db_sample_component["stamp"])
        sample_db["properties"]["stamper"] = db_sample_component["summary"]

        datahandling.save_sample_component_to_file(db_sample_component, sample_component_file)
        datahandling.save_sample_to_file(sample_db, sample_file)
#---Unique to component: end------------------------------------------------------------------------
        open(output, 'w+').close()  # touch file

    except Exception:
        datahandling.write_log(log_out, "Exception in {}\n".format(this_function_name))
        datahandling.write_log(log_err, str(traceback.format_exc()))
        raise Exception
        return 1

    finally:
        datahandling.write_log(log_out, "Done {}\n".format(this_function_name))
        return 0


script__datadump(
    snakemake.output.complete,
    snakemake.params.sample_file,
    snakemake.params.component_file,
    snakemake.params.sample_component_file,
    snakemake.log)

