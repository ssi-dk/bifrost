import functools
import pandas
import os

from bifrostlib import datahandling


# TODO: refactor as requirements_file is the same as component, and component, sample and sample_component are all references to getting the DB entry for each. Currently these references are files and in the future will be ID's
def script__initialization(sample_file, component_file, sample_component_file, output_file, log_out, log_err):
    set_status_to_running(sample_component_file)
    component_db = datahandling.load_component(component_file)
    datahandling.save_component(component_db, component_file)
    if requirements_met(component_file, sample_file, log_out, log_err):
        datahandling.log(log_out, "{}\n{}\n".format(os.getcwd(), output_file))
        with open(str(output_file), "w") as handle:
            handle.write("Requirements met")
    else:
        datahandling.log(log_err, "Requirements not met")
        sample_component_entry = datahandling.load_sample_component(sample_component_file)
        sample_component_entry["status"] = "Requirements not met"
        datahandling.save_sample_component(sample_component_entry, sample_component_file)
    return 0


def set_status_to_running(sample_component_file):
    db_sample_component = datahandling.load_sample_component(sample_component_file)
    db_sample_component["status"] = "Running"
    datahandling.save_sample_component(db_sample_component, sample_component_file)
    return 0


def requirements_met(component_file, sample, log_out, log_err):
    doc_component = datahandling.load_component(component_file)
    db_component_entry = datahandling.load_component(component_id)
    component_file = datahandling.load_yaml(requirements_file)
    sample_name = datahandling.load_yaml(sample)["name"]
    if not passes_check_reads_pipeline(sample, requirements_file, log_err):
        return False
    no_failures = True
    if requirements_file.get('requirements', None) is not None:
        df = pandas.io.json.json_normalize(requirements_file.get('requirements'))  # flattens json into key while maintaining values https://stackoverflow.com/a/41801708
        requirements_dict = df.to_dict(orient='records')[0]

        requirements = []
        for key in requirements_dict:
            values = key.split('.')
            if values[0] == 'components':
                file_location = sample_name + "__" + values[1] + ".yaml"
                requirement = values[2:]  # TODO: add check for no requirements
                expected_value = requirements_dict[key]
                requirements.append([file_location, requirement, expected_value])
            elif values[0] == 'sample':
                file_location = "sample.yaml" 
                requirement = values[1:]  # TODO: add check for no requirements
                expected_value = requirements_dict[key]
                requirements.append([file_location, requirement, expected_value])
            elif values[0] == 'run':
                file_location = "run.yaml"
                requirement = values[1:]  # TODO: add check for no requirements
                expected_value = requirements_dict[key]
                requirements.append([file_location, requirement, expected_value])
            else:
                datahandling.log(log_err, "Improper requirement {}".format(key))

        for requirement in requirements:
            # requirements are in form [file_path, [keys,key2,...], value_of_key (optional otherwise None)]
            file_location = requirement[0]
            keys = requirement[1]
            desired_value = requirement[2]

            db = datahandling.load_yaml(file_location)

            """
            What this does is run dict.get interatively on the db based on the keys till it can 't go
            deeper then returns the value or None if it couldn' t reach that level. While used with dict.get
            any function can be passed
            https://pymotw.com/3/functools/
            """
            actual_value = functools.reduce(dict.get, keys, db)

            """
            Check has been adjusted to check for a list to allow multiple potential options to match
            """
            if not isinstance(desired_value, list):
                desired_value = [desired_value]

            # As it should be a list because of the last line.
            if actual_value is not None:
                # Not sure why desired is [None] instead of None
                if desired_value != [None]:
                    if actual_value in desired_value:
                        datahandling.log(log_err, "Found required entry (value checked) for\ndb: {}\nentry: {}\n".format(":".join(keys), db))
                    else:
                        datahandling.log(log_err, "Requirements not met for\ndb: {}\nentry: {}\ndesired_entry: {}\n".format(":".join(keys), db, desired_value))
                        no_failures = False
                else:
                    datahandling.log(log_err, "Found required entry (value not checked) for\ndb: {}\nentry: {}\n".format(":".join(keys), db))
            else:
                datahandling.log(log_err, "Requirements not met for\ndb: {}\nentry: {}\n".format(file_location, ":".join(keys)))
                no_failures = False
    if no_failures:
        return True
    else:
        return False


def passes_check_reads_pipeline(sample, requirements_file, log_err):
    """
    Checks if the component is a pipeline. In that case it will require reads to be present
    so the component can run.
    """
    sample_db = datahandling.load_yaml(sample)
    if requirements_file["type"] == "pipeline":
        if "reads" not in sample_db:
            datahandling.log(
                log_err, "Pipeline component can't run on a sample with no reads. db:{}".format(sample_db))
            return False
    return True
