import sys
import datetime
from bifrostlib import datahandling

def test__sample__has_reads_files(sampleComponentObj):
    try:
        this_function_name = sys._getframe().f_code.co_name
        summary, results, file_path, key = sampleComponentObj.start_data_extraction()
        test = datahandling.stamperTestObj(this_function_name, "No reads", "core facility")
        if sampleComponentObj.get_reads() == ("", ""):
            test.set_status_and_reason("fail", "Read path is empty")
        else:
            test.set_status_and_reason("pass", "")
    except KeyError as e:
        test.set_status_and_reason("fail", "Database KeyError {} in function {}: ".format(e.args[0], this_function_name))
    finally:
        summary[this_function_name] = test.as_dict()
        return (summary, results)


def test__size_check__has_min_reads(sampleComponentObj):
    try:
        this_function_name = sys._getframe().f_code.co_name
        summary, results, file_path, key = sampleComponentObj.start_data_extraction()
        size_check = sampleComponentObj.get_sample_properties_by_category("size_check")
        test = datahandling.stamperTestObj(this_function_name, "Less than min reads", "core facility")
        test.set_value(size_check["has_min_num_of_reads"])
        if test.get_value() == False:
            test.set_status_and_reason("fail", "Less than min reads")
        else:
            test.set_status_and_reason("pass", "")
    except KeyError as e:
        test.set_status_and_reason("fail", "Database KeyError {} in function {}: ".format(e.args[0], this_function_name))
    finally:
        summary[this_function_name] = test.as_dict()
        return (summary, results)


def test__species_detection__main_species_level(sampleComponentObj):
    try:
        this_function_name = sys._getframe().f_code.co_name
        summary, results, file_path, key = sampleComponentObj.start_data_extraction()
        test = datahandling.stamperTestObj(this_function_name, "Multiple species detected", "supplying lab")
        options = sampleComponentObj.get_options()
        species_detection = sampleComponentObj.get_sample_properties_by_category("species_detection")
        test.set_value(species_detection["percent_classified_species_1"] + species_detection["percent_unclassified"])
        if test.get_value() < options["min_species"]:
            test.set_status_and_reason("fail", "Value ({}) is below threshold ({})".format(test.get_value(), options["min_species"]))
        else:
            test.set_status_and_reason("pass", "")
    except KeyError as e:
        test.set_status_and_reason("fail", "Database KeyError {} in function {}: ".format(e.args[0], this_function_name))
    finally:
        summary[this_function_name] = test.as_dict()
        return (summary, results)


def test__species_detection__unclassified_level(sampleComponentObj):
    try:
        this_function_name = sys._getframe().f_code.co_name
        summary, results, file_path, key = sampleComponentObj.start_data_extraction()
        test = datahandling.stamperTestObj(this_function_name, "High unclassified", "supplying lab")
        options = sampleComponentObj.get_options()
        species_detection = sampleComponentObj.get_sample_properties_by_category("species_detection")
        test.set_value(species_detection["percent_unclassified"])
        if test.get_value() >= options["max_unclassified"]:
            test.set_status_and_reason("fail", "Value ({}) is above threshold ({})".format(test.get_value(), options["max_unclassified"]))
        else:
            test.set_status_and_reason("pass", "")
    except KeyError as e:
        test.set_status_and_reason("fail", "Database KeyError {} in function {}: ".format(e.args[0], this_function_name))
    finally:
        summary[this_function_name] = test.as_dict()
        return (summary, results)


def test__component__species_in_db(sampleComponentObj):
    try:
        this_function_name = sys._getframe().f_code.co_name
        summary, results, file_path, key = sampleComponentObj.start_data_extraction()
        test = datahandling.stamperTestObj(this_function_name, "No species submitted - using default values", "supplying lab")
        options = sampleComponentObj.get_options()
        species_detection = sampleComponentObj.get_sample_properties_by_category("species_detection")
        test.set_value(species_detection["species"])
        if test.get_value() not in options["species_qc_value_mapping"]:
            test.set_status_and_reason("fail", "Detected species not in bifrost db. Can't estimate proper QC values.")
        else:
            test.set_status_and_reason("pass", "")
    except KeyError as e:
        test.set_status_and_reason("fail", "Database KeyError {} in function {}: ".format(e.args[0], this_function_name))
    finally:
        summary[this_function_name] = test.as_dict()
        return (summary, results)


def test__sample__species_provided_is_detected(sampleComponentObj):
    try:
        this_function_name = sys._getframe().f_code.co_name
        summary, results, file_path, key = sampleComponentObj.start_data_extraction()
        test = datahandling.stamperTestObj(this_function_name, "Detected species mismatch", "supplying lab")
        options = sampleComponentObj.get_options()
        species_detection = sampleComponentObj.get_sample_properties_by_category("species_detection")
        sample_info = sampleComponentObj.get_sample_properties_by_category("sample_info")
        test.set_value(sample_info.get("provided_species", None))
        if test.get_value() is None:
            test.set_status_and_reason("pass", "No submitted species")
        elif test.get_value() not in options["species_qc_value_mapping"]:
            test.set_status_and_reason("pass", "Submitted species not in db")
        elif species_detection["species"] != test.get_value():
            test.set_status_and_reason("fail", "Detected species ({}) different than expected ({})".format(species_detection["species"], test.get_value()))
        else:
            test.set_status_and_reason("pass", "")
    except KeyError as e:
        test.set_status_and_reason("fail", "Database KeyError {} in function {}: ".format(e.args[0], this_function_name))
    finally:
        summary[this_function_name] = test.as_dict()
        return (summary, results)


def test__denovo_assembly__genome_size_at_1x(sampleComponentObj):
    try:
        this_function_name = sys._getframe().f_code.co_name
        summary, results, file_path, key = sampleComponentObj.start_data_extraction()
        test = datahandling.stamperTestObj(this_function_name, "Atypical genome size (1x)", "supplying lab")
        options = sampleComponentObj.get_options()
        denovo_assembly = sampleComponentObj.get_sample_properties_by_category("denovo_assembly")
        species_detection = sampleComponentObj.get_sample_properties_by_category("species_detection")
        test.set_value(denovo_assembly["bin_length_at_1x"])
        species = species_detection["species"]
        if species not in options["species_qc_value_mapping"]:
            species = "default"
        if options["species_qc_value_mapping"][species]["min_length"] < test.get_value() < options["species_qc_value_mapping"][species]["max_length"]:
            test.set_status_and_reason("pass", "")
        else:
            test.set_status_and_reason("fail", "Value ({}) below or above expected ({}, {})".format(test.get_value(), options["species_qc_value_mapping"][species]["min_length"], options["species_qc_value_mapping"][species]["max_length"]))
    except KeyError as e:
        test.set_status_and_reason("fail", "Database KeyError {} in function {}: ".format(e.args[0], this_function_name))
    finally:
        summary[this_function_name] = test.as_dict()
        return (summary, results)


def test__denovo_assembly__genome_size_at_10x(sampleComponentObj):
    try:
        this_function_name = sys._getframe().f_code.co_name
        summary, results, file_path, key = sampleComponentObj.start_data_extraction()
        test = datahandling.stamperTestObj(this_function_name, "Atypical genome size (10x)", "supplying lab")
        options = sampleComponentObj.get_options()
        denovo_assembly = sampleComponentObj.get_sample_properties_by_category("denovo_assembly")
        species_detection = sampleComponentObj.get_sample_properties_by_category("species_detection")
        test.set_value(denovo_assembly["bin_length_at_10x"])
        species = species_detection["species"]
        if species not in options["species_qc_value_mapping"]:
            species = "default"
        if options["species_qc_value_mapping"][species]["min_length"] < test.get_value() < options["species_qc_value_mapping"][species]["max_length"]:
            test.set_status_and_reason("pass", "")
        else:
            test.set_status_and_reason("fail", "Value ({}) below or above expected ({}, {})".format(test.get_value(), options["species_qc_value_mapping"][species]["min_length"], options["species_qc_value_mapping"][species]["max_length"]))
    except KeyError as e:
        test.set_status_and_reason("fail", "Database KeyError {} in function {}: ".format(e.args[0], this_function_name))
    finally:
        summary[this_function_name] = test.as_dict()
        return (summary, results)


def test__denovo_assembly__genome_size_difference_1x_10x(sampleComponentObj):
    try:
        this_function_name = sys._getframe().f_code.co_name
        summary, results, file_path, key = sampleComponentObj.start_data_extraction()
        test = datahandling.stamperTestObj(this_function_name, "Atypical genome size difference (1x - 10x)", "supplying lab")
        options = sampleComponentObj.get_options()
        denovo_assembly = sampleComponentObj.get_sample_properties_by_category("denovo_assembly")
        test.set_value(denovo_assembly["bin_contigs_at_1x"] - denovo_assembly["bin_contigs_at_10x"])
        if test.get_value() < options["max_size_difference_for_1x_and_10x"]:
            test.set_status_and_reason("pass", "")
        else:
            test.set_status_and_reason("fail", "Value ({}) above expected ({})".format(test.get_value(), options["max_size_difference_for_1x_and_10x"]))
    except KeyError as e:
        test.set_status_and_reason("fail", "Database KeyError {} in function {}: ".format(e.args[0], this_function_name))
    finally:
        summary[this_function_name] = test.as_dict()
        return (summary, results)


def test__denovo_assembly__genome_average_coverage(sampleComponentObj):
    try:
        this_function_name = sys._getframe().f_code.co_name
        summary, results, file_path, key = sampleComponentObj.start_data_extraction()
        test = datahandling.stamperTestObj(this_function_name, "Atypical genome coverage", "supplying lab")
        options = sampleComponentObj.get_options()
        denovo_assembly = sampleComponentObj.get_sample_properties_by_category("denovo_assembly")
        test.set_value(denovo_assembly["bin_coverage_at_1x"])
        if test.get_value() < options["average_coverage_fail"]:
            test.set_status_and_reason("fail", "Lack of reads ({} < {})".format(test.get_value(), options["average_coverage_fail"]))
        elif test.get_value() < options["average_coverage_low"]:
            test.set_status_and_reason("fail", "Not enough reads ({} < {})".format(test.get_value(), options["average_coverage_low"]))
        elif test.get_value() < options["average_coverage_warn"]:
            test.set_status_and_reason("fail", "Low reads ({} < {})".format(test.get_value(), options["average_coverage_warn"]))
            test.set_effect("supplying lab")
        else:
            test.set_status_and_reason("pass", "")
    except KeyError as e:
        test.set_status_and_reason("fail", "Database KeyError {} in function {}: ".format(e.args[0], this_function_name))
    finally:
        summary[this_function_name] = test.as_dict()
        return (summary, results)


def evaluate_tests_and_stamp(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction()
    species_detection = sampleComponentObj.get_sample_properties_by_category("species_detection")
    sample_info = sampleComponentObj.get_sample_properties_by_category("sample_info")
    core_facility = False
    supplying_lab = False
    for test in summary:
        if summary[test]["status"] == "fail" or summary[test]["status"] == "undefined":
            if summary[test]["effect"] == "supplying lab":
                supplying_lab = True
            elif summary[test]["effect"] == "core facility":
                core_facility = True
    if (sample_info.get("provided_species", None) == species_detection["detected_species"] and
        summary["test__denovo_assembly__genome_average_coverage"]["status"] == "fail" and \
        summary["test__denovo_assembly__genome_average_coverage"]["effect"] == "supplyinh lab" and \
        summary["test__denovo_assembly__genome_size_difference_1x_10x"]["status"] == "fail" and \
        summary["test__sample__species_provided_is_detected"]["status"] == "pass" and \
        summary["test__denovo_assembly__genome_size_at_1x"]["status"] == "pass"):
            core_facility = True
    action = "OK"
    status = "pass"
    if supplying_lab:
        status = "fail"
        action = "supplying lab"
    if core_facility:
        status = "fail"
        action = "core facility"
    summary["stamp"] = {
        "display_name": "ssi_stamper",
        "name": "ssi_stamper",
        "status": status,
        "value": action,
        "date": datetime.datetime.utcnow(),
        "reason": ""
    }
    return (summary, results)


def generate_report(sampleComponentObj):
    summary, results, file_path, key = sampleComponentObj.start_data_extraction()
    data = []
    for test in summary:
        data.append({"test": "{}: {}:{}:{}".format(summary[test]["display_name"],
                                                   summary[test]["status"],
                                                   summary[test]["value"],
                                                   summary[test]["reason"])})
    return data


def datadump(sampleComponentObj, log):
    sampleComponentObj.start_data_dump(log=log)
    sampleComponentObj.run_data_dump_on_function(test__sample__has_reads_files, log=log)
    sampleComponentObj.run_data_dump_on_function(test__size_check__has_min_reads, log=log)
    sampleComponentObj.run_data_dump_on_function(test__species_detection__main_species_level, log=log)
    sampleComponentObj.run_data_dump_on_function(test__species_detection__unclassified_level, log=log)
    sampleComponentObj.run_data_dump_on_function(test__component__species_in_db, log=log)
    sampleComponentObj.run_data_dump_on_function(test__sample__species_provided_is_detected, log=log)
    sampleComponentObj.run_data_dump_on_function(test__denovo_assembly__genome_size_at_1x, log=log)
    sampleComponentObj.run_data_dump_on_function(test__denovo_assembly__genome_size_at_10x, log=log)
    sampleComponentObj.run_data_dump_on_function(test__denovo_assembly__genome_size_difference_1x_10x, log=log)
    sampleComponentObj.run_data_dump_on_function(test__denovo_assembly__genome_average_coverage, log=log)
    sampleComponentObj.run_data_dump_on_function(evaluate_tests_and_stamp, log=log)
    sampleComponentObj.end_data_dump(generate_report_function=generate_report, log=log)


datadump(
    snakemake.params.sampleComponentObj,
    snakemake.log)
