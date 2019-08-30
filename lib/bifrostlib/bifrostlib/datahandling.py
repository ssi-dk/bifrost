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
import sys
import subprocess


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
class stamperTestObj:
    def __init__(self, function_name, display_name, effect, log=None):
        self.name = function_name,
        self.display_name = display_name,
        self.effect = effect
        self.value = ""
        self.status = ""
        self.reason = ""
        self.log = log
        self.write_log_out("Running {}".format(self.name))

    def set_value(self, value):
        if isinstance(value, float):
            value = round(value, 3)
        self.value = value

    def get_value(self):
        return self.value

    def set_status_and_reason(self, status, reason):
        self.status = status
        self.reason = self.reason

    def write_log_out(self, content):
        if self.log is not None:
            with open(self.log.log_out, "a+") as file_handle:
                file_handle.write(content)
        else:
            sys.stdout.write(content)

class SampleComponentObj:
    def __init__(self):
        self.sample_id = None
        self.component_id = None

    def load(self, sample_id, component_id):
        self.sample_id = sample_id
        self.component_id = component_id
        self.sample_db = get_sample(sample_id=self.sample_id)
        self.component_db = get_component(component_id=self.component_id)
        self.sample_component_db = get_sample_component(sample_id=self.sample_id, component_id=self.component_id)
        self.sample_component_id = self.sample_component_db["_id"]
        return (self.sample_db["name"], self.component_db["name"], self.component_db["dockerfile"], self.component_db["options"], self.component_db["resources"])

    def start_data_extraction(self, file_location=None):
        summary = self.sample_component_db["properties"]["summary"]
        results = self.sample_component_db["results"]
        file_path = None
        key = None
        if file_location is not None:
            file_path = os.path.join(self.get_component_name(), file_location)
            key = self.get_file_location_key(file_location)
            results[key] = {}
        return summary, results, file_path, key

    def get_file_location_key(self, file_location):
        file_path = os.path.join(self.get_component_name(), file_location)
        key = file_path.replace(".", "_").replace("$", ".")
        return key

    def get_sample_properties_by_category(self, category):
        if category in self.sample_db["properties"]:
            return self.sample_db["properties"][category].get("summary", None)
        else:
            return None

    def get_reads(self):
        datafiles = self.get_sample_properties_by_category("datafiles")
        if "paired_reads" in datafiles:
            return (datafiles["paired_reads"][0], datafiles["paired_reads"][1])
        else:
            return ("/dev/null", "/dev/null")

    def get_resources(self):
        return self.component_db["resources"]

    def get_options(self):
        return self.component_db["options"]

    def check_requirements(self, output_file="requirements_met", log=None):
        no_failures = True
        if self.component_db["requirements"] is not None:
            requirements = pandas.io.json.json_normalize(self.component_db["requirements"], sep=".").to_dict(orient='records')[0]  # a little loaded of a line, get requirements from component_db, use the pandas json function to turn it into a 2d dataframe, then convert that to a dict of known depth 2, 0 is for our 1 and only sheet
            for requirement in requirements:
                category = requirement.split(".")[0]
                if category == "sample":
                    field = requirement.split(".")[1:]
                    expected_value = requirements[requirement]
                    if not self.requirement_met(self.sample_db, field, expected_value, log):
                        no_failures = False
                elif category == "component":
                    field = requirement.split(".")[2:]
                    expected_value = requirements[requirement]
                    if not self.requirement_met(self.component_db, field, expected_value, log):
                        no_failures = False
                else:
                    no_failures = False
        if no_failures:
            open(os.path.join(self.component_db["name"], output_file), "w+").close()
        else:
            self.requirements_not_met()

    def requirements_not_met(self):
        sys.stdout.write("Workflow stopped due to requirements\n")
        update_status_in_sample_and_sample_component(self.sample_component_id, "Requirements not met")

    def started(self):
        sys.stdout.write("Workflow processing\n")
        update_status_in_sample_and_sample_component(self.sample_component_id, "Running")

    def success(self):
        sys.stdout.write("Workflow complete\n")
        update_status_in_sample_and_sample_component(self.sample_component_id, "Success")

    def failure(self):
        sys.stdout.write("Workflow error\n")
        update_status_in_sample_and_sample_component(self.sample_component_id, "Failure")

    def start_rule(self, rule_name, log=None):
        self.write_log_out(log, "{} has started\n".format(rule_name))
        return (self.component_db["name"], self.component_db["options"], self.component_db["resources"])

    def rule_run_cmd(self, command, log):
        self.write_log_out(log, "Running:{}\n".format(command))
        command_log_out, command_log_err = subprocess.Popen(command, shell=True).communicate()
        self.write_log_out(log, command_log_out)
        self.write_log_err(log, command_log_err)

    def end_rule(self, rule_name, log=None):
        self.write_log_out(log, "{} has finished\n".format(rule_name))

    def start_data_dump(self, log=None):
        self.sample_component_db["properties"] = {
            "summary": {},
            "component": {
                "_id": self.component_db["_id"]
            }
        }
        self.sample_component_db["results"] = {}
        self.sample_component_db["report"] = self.component_db["db_values_changes"]["sample"]["report"][self.component_db["details"]["category"]]
        self.write_log_out(log, "Starting datadump\n")
        self.save_files_to_sample_component(log)

    def save_files_to_sample_component(self, log=None):
        save_files_to_db(self.component_db["db_values_changes"]["files"], sample_component_id=self.sample_component_db["_id"])
        self.write_log_out(log, "Files saved: {}\n".format(",".join(self.component_db["db_values_changes"]["files"])))

    def get_component_name(self):
        return self.component_db["name"]

    def run_data_dump_on_function(self, data_extraction_function, log=None):
        (self.sample_component_db["properties"]["summary"], self.sample_component_db["results"]) = data_extraction_function(self)

    def end_data_dump(self, output_file="datadump_complete", generate_report_function=lambda x: None, log=None):
        self.sample_db["properties"][self.component_db["details"]["category"]] = self.sample_component_db["properties"]
        report_data = generate_report_function(self)
        if report_data is not None:
            self.sample_db["report"][self.component_db["details"]["category"]]["data"] = report_data
            assert(type(self.sample_db["report"][self.component_db["details"]["category"]]["data"])==list)
        self.write_log_err(log, str(traceback.format_exc()))
        save_sample(self.sample_db)
        self.write_log_out(log, "sample {} saved\n".format(self.sample_db["_id"]))
        save_sample_component(self.sample_component_db)
        self.write_log_out(log, "sample_component {} saved\n".format(self.sample_component_db["_id"]))
        open(os.path.join(self.component_db["name"], output_file), "w+").close()
        self.write_log_out(log, "Done datadump\n")

    def requirement_met(self, db, field, expected_value, log):
        try:
            actual_value = functools.reduce(dict.get, field, db)
            if expected_value is None:
                self.write_log_err(log, "Found required entry (value not checked) for entry: {}\n".format(":".join(field)))
                return True
            elif type(expected_value) is not list:
                expected_value = [expected_value]
            if actual_value in expected_value:
                    self.write_log_err(log, "Found required entry for entry: {} value:{}\n".format(":".join(field), actual_value))
                    return True
            else:
                self.write_log_err(log, "Requirements not met for entry: {} allowed_values: {} value:{}\n".format(":".join(field), expected_value, actual_value))
                return False
        except Exception:
            self.write_log_err(log, "Requirements not met for entry: {}\ndb was:{}\n".format(":".join(field), db))
            self.write_log_err(log, str(traceback.format_exc()))
            return False

    def write_log_out(self, log, content):
        if log is not None:
            self.write_log(log.out_file, content)
        else:
            sys.stdout.write(content)

    def write_log_err(self, log, content):
        if log is not None:
            self.write_log(log.err_file, content)
        else:
            sys.stderr.write(content)

    def write_log(self, log_file, content):
        with open(log_file, "a+") as file_handle:
            file_handle.write(content)


def write_log(log_file, content):
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
    if len(components) == 1:
        return components[0]
    else:
        return temp

def save_sample(sample_db):
    #TODO: This writing to a file is temporary and should be removed with bifrost.smk refactoring
    if "path" in sample_db:
        with open(os.path.join(sample_db["path"], "sample.yaml"), "w") as file_handle:
            yaml.dump(sample_db, file_handle)

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
    #TODO: This writing to a file is temporary and should be removed with bifrost.smk refactoring
    sample_db = get_sample(sample_component_db["sample"]["_id"])
    if "path" in sample_db:
        with open(os.path.join(sample_db["path"], "{}__{}.yaml".format(sample_component_db["sample"]["name"], sample_component_db["component"]["name"])), "w") as file_handle:
            yaml.dump(sample_component_db, file_handle)

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
        if component["_id"] == component_db["_id"]:
            component["status"] = status
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


def load_file_from_db(file_id, save_to_path=None, subpath=False):
    return mongo_interface.load_file_from_db(file_id, save_to_path, subpath)


def recreate_s_c_files(sample_component_id, save_to_path):
    files = mongo_interface.find_files(ObjectId(sample_component_id))
    for f in files:
        load_file_from_db(f._id, subpath=True)

def recreate_yaml(collection, oid, path):
    if collection == "sample_components":
        obj = get_sample_component(oid)
    elif collection == "samples":
        obj = get_sample(oid)
    else:
        raise ValueError("invalid collection")
    with open(path, "w") as file_handle:
        yaml.dump(obj, file_handle)


def test():
    print("Hello")
    return (u'This is the bifrost lib')
