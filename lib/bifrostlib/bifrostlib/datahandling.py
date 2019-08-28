import os
import datetime
import ruamel.yaml
import pandas
import functools
from bson.objectid import ObjectId
from bson.int64 import Int64
from bifrostlib import mongo_interface
import pymongo
import traceback


ObjectId.yaml_tag = u'!bson.objectid.ObjectId'
ObjectId.to_yaml = classmethod(
    lambda cls, representer, node: representer.represent_scalar(
        cls.yaml_tag, u'{}'.format(node)))
ObjectId.from_yaml = classmethod(
    lambda cls, constructor, node: cls(node.value))

Int64.yaml_tag = u'!bson.int64.Int64'
Int64.to_yaml = classmethod(
    lambda cls, representer, node: representer.represent_scalar(
        cls.yaml_tag, u'{}'.format(node)))
Int64.from_yaml = classmethod(
    lambda cls, constructor, node: cls(node.value))

yaml = ruamel.yaml.YAML(typ="safe")
yaml.default_flow_style = False
yaml.register_class(ObjectId)
yaml.register_class(Int64)

"""
Class to be used as a template for rules which require python scripts. Can be tightened up to only
allow access to parts of the document if needed in the future, ie options.
"""
class SampleComponentObj:
    def __init__(self, sample_id, component_id):
        self.sample_id = sample_id
        self.component_id = component_id
        self.load()
    def load(self):
        self.db_sample = get_sample(sample_id=self.sample_id)
        self.db_component = get_component(component_id=self.component_id)
        self.db_sample_component = get_sample_component(sample_id=self.sample_id, component_id=self.component_id)
        self.sample_component_id = self.db_sample_component["_id"]
        self.started()
    def check_requirements(self, log=None, output_file="requirements_met"):
        no_failures = True
        if self.db_component["requirements"] is not None:
            requirements = pandas.io.json.json_normalize(self.db_component["requirements"], sep=".").to_dict(orient='records')[0]  # a little loaded of a line, get requirements from db_component, use the pandas json function to turn it into a 2d dataframe, then convert that to a dict of known depth 2, 0 is for our 1 and only sheet
            for requirement in requirements:
                category = requirement.split(".")[0]
                if category == "sample":
                    field = requirement.split(".")[1:]
                    expected_value = requirements[requirement]
                    if not requirement_met(self.db_sample, field, expected_value, log):
                        no_failures = False
                elif category == "component":
                    field = requirement.split(".")[2:]
                    expected_value = requirements[requirement]
                    if not requirement_met(self.db_component, field, expected_value, log):
                        no_failures = False
                else:
                    write_log_err(log, "Improper requirement {}".format(requirement))
                    no_failures = False
        if no_failures:
            open(os.path.join(self.db_component["name"], output_file), "w+").close()
        else:
            self.requirements_not_met()
    def requirements_not_met(self):
        update_status_in_sample_and_sample_component(self.sample_component_id, "Requirements not met")
    def started(self):
        update_status_in_sample_and_sample_component(self.sample_component_id, "Running")
    def success(self):
        update_status_in_sample_and_sample_component(self.sample_component_id, "Success")
    def failure(self):
        update_status_in_sample_and_sample_component(self.sample_component_id, "Failure")
    def start_rule(self, rule_name, log):
        write_log_out(log, "{} has started\n".format(rule_name))
        return (self.db_sample, self.db_component)
    def end_rule(self, rule_name, log):
        write_log_out(log, "{} has finished\n".format(rule_name))
        return 0
    def start_data_dump(self, log):
        self.db_sample_component["properties"] = {
            "summary": {},
            "component": {
                "_id": self.db_component["_id"]
            }
        }
        self.db_sample_component["results"] = {}
        self.db_sample_component["report"] = self.db_component["db_values_changes"]["sample"]["report"][self.db_component["details"]["category"]]
        write_log_out(log, "Starting datadump\n")
        self.save_files_to_sample_component()
    def save_files_to_sample_component(self):
        try:
            save_files_to_db(self.db_component["db_values_changes"]["files"], sample_component_id=self.db_sample_component["_id"])
            self.write_log_out("Files saved")
        except:
            self.write_log_err(str(traceback.format_exc()))
    def get_summary_and_results(self):
        return (self.db_sample_component["properties"]["summary"], self.db_sample_component["results"])
    def get_component_name(self):
        return self.db_component["name"]
    def run_data_dump_on_function(self, data_extraction_function, log):
        try:
            (self.db_sample_component["properties"]["summary"], self.db_sample_component["results"]) = data_extraction_function(self)
        except Exception:
            write_log_err(log, str(traceback.format_exc()))

    def end_data_dump(self, log, output_file="datadump_complete", generate_report_function=lambda x: None):
        try:
            self.db_sample["properties"][self.db_component["details"]["category"]] = self.db_sample_component["properties"]
            self.db_sample["report"][self.db_component["details"]["category"]] = generate_report_function(self)
            write_log_err(log, str(traceback.format_exc()))
            save_sample_to_file(self.db_sample, self.sample_file)
            write_log_out(log, "sample {} saved\n".format(self.db_sample["_id"]))
            save_sample_component_to_file(self.db_sample_component, self.sample_component_file)
            write_log_out(log, "sample_component {} saved\n".format(self.db_sample_component["_id"]))
            open(os.path.join(self.db_component["name"], output_file), "w+").close()
            write_log_out(log, "Done datadump\n")
        except Exception:
            write_log_err(log, str(traceback.format_exc()))

def requirement_met(db, field, expected_value, log):
    try:
        actual_value = functools.reduce(dict.get, field, db)
        if expected_value is None:
            write_log_err(log, "Found required entry (value not checked) for\ndb: {}\nentry: {}\n".format(":".join(field), db))
            return True
        elif type(expected_value) is not list:
            expected_value = [expected_value]
        if actual_value in expected_value:
                write_log_err(log, "Found required entry (value checked) for\ndb: {}\nentry: {}\n".format(":".join(field), db))
                return True
        else:
            write_log_err(log, "Requirements not met for\ndb: {}\nentry: {}\ndesired_entry: {}\n".format(":".join(field), db, expected_value))
            return False
    except Exception:
        write_log_err(log, "Requirements not met for\ndb: {}\nentry: {}\n".format(db, ":".join(field)))
        write_log_err(log, "Error: " + str(traceback.format_exc()))
        return False
def write_log_out(log, content):
    log(log.out_file, content)
def write_log_err(log, content):
    log(log.err_file, content)





class SampleComponentRuleObj:
    def __init__(self, input, output, sample_file, component_file, log, function_name):
        self.input = input
        self.output = output
        self.sample_file = sample_file
        self.component_file = component_file
        self.log = log
        self.function_name = function_name
        self.load()

    def load(self):
        self.db_sample = load_sample(self.sample_file)
        self.db_component = load_component(self.component_file)
        self.write_log_out("{} has started\n".format(self.function_name))

    def get_sample_and_component_dbs(self):
        return (self.db_sample, self.db_component)

    def get_folder(self):
        return self.db_component["name"]

    def rule_done(self):
        self.write_log_out("{} has finished\n".format(self.function_name))
        return 0

class DatadumpSampleComponentObj:
    def __init__(self, sample_file, component_file, sample_component_file, log):
        self.output_file = "component_complete"
        self.sample_file = sample_file
        self.component_file = component_file
        self.sample_component_file = sample_component_file
        self.log = log
        self.load()

    def load(self):
        self.db_sample = load_sample(self.sample_file)
        self.db_component = load_component(self.component_file)
        self.db_sample_component = load_sample_component(self.sample_component_file)



    def save_files_to_sample_component(self):
        try:
            save_files_to_db(self.db_component["db_values_changes"]["files"], sample_component_id=self.db_sample_component["_id"])
            self.write_log_out("Files saved")
        except:
            self.write_log_err(str(traceback.format_exc()))

    def get_summary_and_results(self):
        return (self.db_sample_component["properties"]["summary"], self.db_sample_component["results"])

    def get_component_name(self):
        return self.db_component["name"]

    def set_summary_and_results_from_function(self, data_extraction_function):
        try:
            (self.db_sample_component["properties"]["summary"], self.db_sample_component["results"]) = data_extraction_function(self)
        except Exception:
            self.write_log_err(str(traceback.format_exc()))

    def save(self, generate_report_function=lambda x: None):
        try:
            self.db_sample["properties"][self.db_component["details"]["category"]] = self.db_sample_component["properties"]
            self.db_sample["report"][self.db_component["details"]["category"]] = generate_report_function(self)
            self.write_log_err(str(traceback.format_exc()))
            save_sample_to_file(self.db_sample, self.sample_file)
            self.write_log_out("sample {} saved\n".format(self.db_sample["_id"]))
            save_sample_component_to_file(self.db_sample_component, self.sample_component_file)
            self.write_log_out("sample_component {} saved\n".format(self.db_sample_component["_id"]))
            open(self.output_path, "w+").close()
            self.write_log_out("Done datadump\n")
        except Exception:
            self.write_log_err(str(traceback.format_exc()))

    def write_log_out(self, content):
        log(self.log.out_file, content)

    def write_log_err(self, content):
        log(self.log.err_file, content)


def log(log_file, content):
    with open(log_file, "a+") as file_handle:
        file_handle.write(content)


def load_config():
    config_file = "config.yaml"
    if not os.path.isfile("config.yaml"):
        config_file = "../config.yaml"
    with open(config_file, "r") as file_handle:
        return yaml.load(file_handle)


def save_run(run_dict, file_yaml):
    config = load_config()
    if config["use_mongodb"]:
        run_dict = mongo_interface.dump_run_info(run_dict)
    with open(file_yaml, "w") as file_handle:
        yaml.dump(run_dict, file_handle)


def load_run(file_yaml):
    if not os.path.isfile(file_yaml):
        return {}
    with open(file_yaml, "r") as file_handle:
        return yaml.load(file_handle)


def save_component(component_dict, file_yaml):
    config = load_config()
    if config["use_mongodb"]:
        component_dict = mongo_interface.dump_component_info(component_dict)
    with open(file_yaml, "w") as file_handle:
        yaml.dump(component_dict, file_handle)


def load_component(file_yaml):
    with open(file_yaml, "r") as file_handle:
        temp = yaml.load(file_handle)
    components = get_components(component_names=[temp["name"]], component_versions=[temp["version"]])
    print(len(components), temp)
    if len(components) == 1:
        return components[0]
    else:
        return temp

def save_sample(sample_db):
    return mongo_interface.dump_sample_info(sample_db)

def save_sample_to_file(sample_dict, file_yaml):
    config = load_config()
    if (len(file_yaml.split("/")) >= 2 and
            os.path.isdir("/".join(file_yaml.split("/")[:-1])) is False):
        os.mkdir("/".join(file_yaml.split("/")[:-1]))
    if config["use_mongodb"]:
        sample_dict = mongo_interface.dump_sample_info(sample_dict)
    with open(file_yaml, "w") as file_handle:
        yaml.dump(sample_dict, file_handle)


def load_sample(file_yaml):
    if not os.path.isfile(file_yaml):
        return {}
    with open(file_yaml, "r") as file_handle:
        content = yaml.load(file_handle)
        if content is None:
            return {}
        else:
            return content

def save_sample_component(sample_component_db):
    return mongo_interface.dump_sample_component_info(sample_component_db)

def save_sample_component_to_file(sample_component_dict, file_yaml):
    config = load_config()
    if (len(file_yaml.split("/")) >= 2 and
            os.path.isdir("/".join(file_yaml.split("/")[:-1])) is False):
        os.mkdir("/".join(file_yaml.split("/")[:-1]))
    if config["use_mongodb"]:
        sample_component_dict = mongo_interface.dump_sample_component_info(
            sample_component_dict)
    with open(file_yaml, "w") as file_handle:
        yaml.dump(sample_component_dict, file_handle)


def load_sample_component(file_yaml):
    if not os.path.isfile(file_yaml):
        return {}
    with open(file_yaml, "r") as file_handle:
        content = yaml.load(file_handle)
        if content is None:
            return {}
        else:
            return content


def sample_component_success(file_yaml):
    sample_component_dict = load_sample_component(file_yaml)
    if sample_component_dict["status"] != "Success":
        return False
    else:
        return True


def update_status_in_sample_and_sample_component(sample_component_id, status):
    sample_component_db = get_sample_component(sample_component_id=sample_component_id)
    component_db = get_component(component_id=sample_component_db["component"]["_id"])
    sample_db = get_sample(sample_id=sample_component_db["sample"]["_id"])
    sample_component_db["status"] = status
    status_set = False
    for component in sample_db["components"]:
        if sample_db["components"][component]["_id"] == component_db["_id"]:
            sample_db["components"][component]["status"] = status
            status_set = True
    if not status_set:
        sample_db["components"].append([{"_id":component_db["_id"], "name":component_db["name"], "status":status}])
    save_sample(sample_db)
    save_sample_component(sample_component_db)


def update_sample_component_failure(file_yaml):
    sample_component_dict = load_sample_component(file_yaml)
    if sample_component_dict["status"] != "Requirements not met":
        sample_component_dict["status"] = "Failure"
    save_sample_component_to_file(sample_component_dict, file_yaml)


def save_yaml(content_dict, file_yaml):
    with open(file_yaml, "w") as output:
        yaml.dump(content_dict, output)


def load_yaml(file_yaml):
    if not os.path.isfile(file_yaml):
        return {}
    with open(file_yaml, "r") as file_handle:
        return yaml.load(file_handle)


def read_buffer(file_path):
    with open(file_path, "r") as file_handle:
        buffer = file_handle.read()
    return buffer


def datadump_template(extraction_callback, db, temp_data=None, key=None, file_path=None):
    try:
        if file_path is not None:
            if os.path.isfile(file_path):
                if key is None:
                    key = file_path.replace(".", "_").replace("$", "_")  # $ and . are special characters to be avoided in keys
        if key is not None:
            db["results"][key] = {}
        db = extraction_callback(db, file_path=file_path, key=key, temp_data=temp_data)

    except Exception:
        print(traceback.format_exc())
        if key is not None:
            db["results"][key]["status"] = "datadumper error"
        raise Exception

    return db


# /runs
def post_run(run):
    # Used only by test suite for now.
    # NOTE: dump_run_info acts like a PUT
    return mongo_interface.dump_run_info(run)


def get_runs(run_id=None,
             names=None, sample_id=None):
    if run_id is not None:
        run_id = ObjectId(run_id)
    if sample_id is not None:
        sample_id = ObjectId(sample_id)
    return mongo_interface.get_runs(names=names,
                                    sample_id=sample_id,
                                    run_id=run_id)


def delete_run(name=None, run_id=None):
    if run_id is not None:
        run_db = get_runs(run_id=run_id)[0]
    elif name is not None:
        run_db = get_runs(names=[name])[0]

    for sample in run_db["samples"]:
        sample_runs = get_runs(sample_id=str(sample["_id"]))
        if len(sample_runs) == 1:
            delete_sample(str(sample["_id"]))
    
    for component in run_db["components"]:
        samples_with_components = get_samples(
            component_ids=[component["_id"]])
        if len(samples_with_components) == 0:
            delete_component(str(component["_id"]))
    
    deleted = mongo_interface.delete_run(run_db["_id"])
    return deleted


def get_run_export(names=None):
    """
    Export runs
    """
    run_dicts = {}
    for run_name in names:
        run_db = get_runs(names=[run_name])[0]
        component_ids = set()
        if "components" in run_db:
            for comp in run_db["components"]:
                component_ids.add(comp["_id"])
        sample_ids = [str(s["_id"]) for s in run_db["samples"]]
        samples = get_samples(sample_ids=sample_ids)
        for sample in samples:
            if "components" in sample:
                for comp in sample["components"]:
                    component_ids.add(comp["_id"])
        components = get_components(component_ids=list(component_ids))
        sample_components = get_sample_components(sample_ids=sample_ids)
        run_dicts[run_name] = {
            "components": components,
            "samples": samples,
            "sample_components": sample_components,
            "runs": [run_db]
        }
    return run_dicts


def post_run_export(import_dict):
    """
    Import runs
    """
    imported = 0

    for import_run in import_dict.values():

        for c in import_run["components"]:
            try:
                get_c = post_component(c)
                if get_c is not None:
                    imported += 1
            except pymongo.errors.PyMongoError as e:
                print("Error importing component with id {}.\n{}".format(c["_id"], c))
                print("Error: \n{}\n\n".format(str(e)))
            

        for s in import_run["samples"]:
            try:
                get_s = post_sample(s)
                if get_s is not None:
                    imported += 1
            except pymongo.errors.PyMongoError as e:
                print("Error importing sample with id {}.\n{}".format(
                    s["_id"], s))
                print("Error: \n{}\n\n".format(str(e)))
            

        for s_c in import_run["sample_components"]:
            try:
                get_s_c = post_sample_component(s_c)
                if get_s_c is not None:
                    imported += 1
            except pymongo.errors.PyMongoError as e:
                print("Error importing sample_component with id {}.\n{}".format(
                    s_c["_id"], s_c))
                print("Error: \n{}\n\n".format(str(e)))

        for r in import_run["runs"]:
            try:
                get_r = post_run(r)
                if get_r is not None:
                    imported += 1
            except pymongo.errors.PyMongoError as e:
                print("Error importing run with id {}.\n{}".format(
                    r["_id"], r))
                print("Error: \n{}\n\n".format(str(e)))

        return imported

# /samples
def get_samples(sample_ids=None, run_names=None, component_ids=None):
    if sample_ids is not None:
        sample_ids = [ObjectId(id) for id in sample_ids]
    if component_ids is not None:
        component_ids = [ObjectId(id) for id in component_ids]
    return mongo_interface.get_samples(sample_ids=sample_ids,
                                       run_names=run_names,
                                       component_ids=component_ids)


def get_sample(sample_id=None):
    if sample_id is not None:
        sample_id = [ObjectId(sample_id)]
    return next(iter(mongo_interface.get_samples(sample_ids=sample_id)), None)


def post_sample(sample):
    return mongo_interface.dump_sample_info(sample)


def delete_sample(sample_id):
    """
    Delete sample, sample_components and sample from runs.
    Only from database.
    """
    delete_sample_component(sample_id=sample_id)
    delete_sample_from_runs(sample_id=sample_id)
    deleted = mongo_interface.delete_sample(sample_id=ObjectId(sample_id))
    return deleted


#  /sample/{id}/sample_components


def get_sample_components(sample_component_ids=None, sample_ids=None, component_ids=None, component_names=None):
    # Should be smarter
    if sample_component_ids is not None:
        sample_component_ids = [ObjectId(id) for id in sample_component_ids]
    if sample_ids is not None:
        sample_ids = [ObjectId(id) for id in sample_ids]
    if component_ids is not None:
        component_ids = [ObjectId(id) for id in component_ids]
    return mongo_interface.get_sample_components(
        sample_component_ids=sample_component_ids,
        sample_ids=sample_ids,
        component_names=component_names)


def get_sample_component(sample_component_id=None, sample_id=None, component_id=None):
    if sample_component_id is not None:
        sample_component_id = [ObjectId(sample_component_id)]
    if sample_id is not None:
        sample_id = [ObjectId(sample_id)]
    if component_id is not None:
        component_id = [ObjectId(component_id)]
    return next(iter(mongo_interface.get_sample_components(sample_component_ids=sample_component_id, sample_ids=sample_id, component_ids=component_id)), None)


def post_sample_component(sample_component):
    return mongo_interface.dump_sample_component_info(sample_component)


def delete_sample_component(s_c_id=None, sample_id=None):
    if s_c_id is not None:
        s_c_id = ObjectId(s_c_id)
    if sample_id is not None:
        sample_id = ObjectId(sample_id)
    return mongo_interface.delete_sample_component(s_c_id, sample_id)


# /component

def get_components(component_ids=None, component_names=None, component_versions=None):
    """
    Get components from db
    """
    if component_ids is not None:
        component_ids = list(map(ObjectId, component_ids))
    return mongo_interface.get_components(component_ids=component_ids, component_names=component_names, component_versions=component_versions)


def get_component(component_id=None):
    if component_id is not None:
        component_id = [ObjectId(component_id)]
    return next(iter(mongo_interface.get_components(component_ids=component_id)), None)


def post_component(component):
    return mongo_interface.dump_component_info(component)


def delete_component(component_id):
    return mongo_interface.delete_component(ObjectId(component_id))


# species
def get_mlst_species_DB(file_yaml):
    with open(file_yaml, "r") as file_handle:
        sample = yaml.load(file_handle)
        species = sample["properties"].get("species", None)
        return mongo_interface.query_mlst_species(species)


def load_species(ncbi_species):
    return mongo_interface.query_species(ncbi_species)


def get_ncbi_species(species_entry):
    return mongo_interface.query_ncbi_species(species_entry)

# Refactor this


def load_samples_from_runs(run_ids=None, names=None):
    if run_ids is not None:
        if type(run_ids) == str:
            run_ids = [run_ids]
        run_ids = [ObjectId(id) for id in run_ids]
    res = mongo_interface.load_samples_from_runs(run_ids, names)
    sample_ids = set()  # avoid dupes
    for run in res:
        sample_ids.update(list(map(lambda x: x["_id"], run["samples"])))
    return list(sample_ids)


def delete_sample_from_runs(sample_id=None):
    if sample_id is not None:
        sample_id = ObjectId(sample_id)
    return mongo_interface.delete_sample_from_runs(sample_id)



def save_files_to_db(file_paths, sample_component_id):
    if file_paths is None:
        file_paths = []
    file_ids = []
    for file_path in file_paths:
        file_ids.append(mongo_interface.save_file_to_db(
            sample_component_id, file_path))
    return file_ids


def load_file_from_db(file_id, save_to_path=None):
    return mongo_interface.load_file_from_db(file_id, save_to_path)

def test():
    print("Hello")
    return (u'This is the bifrost lib')
