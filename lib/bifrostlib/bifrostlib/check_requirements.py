import functools
import pandas
import os
import traceback

from bifrostlib import datahandling

def script__initialization(sample_file, component_file, sample_component_file, output_file, log_out, log_err):
    set_status_to_running(sample_component_file)
    component_db = datahandling.load_component(component_file)
    if all_requirements_met(component_file, sample_file, log_out, log_err):
        datahandling.log(log_out, "{}\n{}\n".format(os.getcwd(), output_file))
        with open(str(output_file), "w") as handle:
            handle.write("Requirements met")
    else:
        datahandling.log(log_out, "Requirements not met")
        sample_component_entry = datahandling.load_sample_component(sample_component_file)
        sample_component_entry["status"] = "Requirements not met"
        datahandling.save_sample_component(sample_component_entry, sample_component_file)
    return 0


def set_status_to_running(sample_component_file):
    db_sample_component = datahandling.load_sample_component(sample_component_file)
    db_sample_component["status"] = "Running"
    datahandling.save_sample_component(db_sample_component, sample_component_file)
    return 0


def requirement_met(db, field, expected_value, log_out, log_err):
    try:
        actual_value = functools.reduce(dict.get, field, db)
        if expected_value is None:
            datahandling.log(log_err, "Found required entry (value not checked) for\ndb: {}\nentry: {}\n".format(":".join(field), db))
            return True
        elif type(expected_value) is list:
            if actual_value in expected_value:
                datahandling.log(log_err, "Found required entry (value checked) for\ndb: {}\nentry: {}\n".format(":".join(field), db))
                return True
            else:
                datahandling.log(log_err, "Requirements not met for\ndb: {}\nentry: {}\ndesired_entry: {}\n".format(":".join(field), db, expected_value))
                return False
    except Exception:
        datahandling.log(log_err, "Requirements not met for\ndb: {}\nentry: {}\n".format(db, ":".join(field)))
        datahandling.log(log_err, "Error: " + str(traceback.format_exc()))
        return False 


def all_requirements_met(component_file, sample_file, log_out, log_err):
    db_component = datahandling.load_component(component_file)
    db_sample = datahandling.load_sample(sample_file)

    no_failures = True
    db_component["requirements"] = db_component.get("requirements", None)
    if db_component["requirements"] is not None:
        requirements = pandas.io.json.json_normalize(db_component["requirements"], sep=".").to_dict(orient='records')[0]  # a little loaded of a line, get requirements from db_component, use the pandas json function to turn it into a 2d dataframe, then convert that to a dict of known depth 2, 0 is for our 1 and only sheet
        for requirement in requirements:
            category = requirement.split(".")[0]
            if category == "sample":
                field = requirement.split(".")[1:]
                expected_value = requirements[requirement]
                if not requirement_met(db_sample, field, expected_value, log_out, log_err):
                    no_failures = False
            elif category == "component":
                sample_component_file = db_sample["name"] + "__" + requirement.split(".")[1] + ".yaml"
                db_sample_component = datahandling.load_sample_component(sample_component_file)
                field = requirement.split(".")[2:]
                expected_value = requirements[requirement]
                if not requirement_met(db_sample_component, field, expected_value, log_out, log_err):
                    no_failures = False
            else:
                datahandling.log(log_err, "Improper requirement {}".format(requirement))
                no_failures = False
    else:
        datahandling.log(log_err, "No requirements {}".format(requirement))
    if no_failures:
        return True
    else:
        return False
