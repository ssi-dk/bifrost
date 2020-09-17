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
from typing import List, Set, Dict, Tuple, Optional


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
        self.name = function_name
        self.display_name = display_name
        self.effect = effect
        self.value = ""
        self.status = ""
        self.reason = ""
        self.log = log
        self.write_log_out("Running {}\n".format(function_name))

    def set_value(self, value):
        if isinstance(value, float):
            value = round(value, 3)
        self.value = value

    def get_value(self):
        return self.value

    def set_effect(self, effect):
        self.effect = effect

    def set_status_and_reason(self, status, reason):
        self.status = status
        self.reason = reason

    def as_dict(self):
        test = {
            "name": self.name,
            "display_name": self.display_name,
            "effect": self.effect,
            "value": self.value,
            "status": self.status,
            "reason": self.reason,
        }
        return test

    def write_log_out(self, content):
        if self.log is not None:
            with open(self.log.log_out, "a+") as file_handle:
                file_handle.write(content)
        else:
            sys.stdout.write(content)


# class Component():
#     def __init__(self) -> None:
#         self._dict = None
#         if _id != None:
#             components = get_components(_id)
#             if len(components) == 1:
#                 self._dict = components[0]
#         elif:

class Category:
    """
    Category Object for use in Sample Object and Sample Component Objects
    A Category is a resulting value from a component which other components may also do
    an example of this would be 2 different denovo assembly pipelines, both provide contigs but
    in different ways. They would both be the same category so the user could pull common values
    with the caveat of different pipelines are not direct apple to apple comparisons.
        - Summary, should be the same summary as created by the component for the sample_component
        - Component, a reference id to the component which created the category.
    """
    def __init__(self, name: str) -> None:
        self.name = name # Check for name in list?
        self._dict = None
        self._dict = {
            "summary": {},
        }

    def get_name(self) -> str:
        """Returns name of object"""
        return self.name 

    def get(self, key: str) -> object:
        """Returns value of the associated key in object (dict)"""
        return self._dict[key]

    def set_summary(self, summary: dict) -> None:
        """Sets summary value in object (dict)"""
        self._dict["summary"] = summary

    def set_component(self, component_id: str) -> None:
        """Sets component value in object (dict)"""
        self._dict["component"] = {"_id": component_id}

    def display(self) -> dict:
        """Return a copy of the object (dict)"""
        return self._dict.copy()

class Sample:
    """
    Sample Object for saving in bifrost DB document "samples"
    A Sample is a base unit in the bifrost DB. It contains information regarding:
        - Data files
        - Sample info (light metadata)
        - Components that it the sample has been run against and status
        - Report, or display information for web reporter
        - Shared category properties for components with similar function
        - Metadata (for the DB)
    When Sample.save() is ran the values are committed to the DB. 
    """
    def __init__(self, _id: str = None, name: str = None) -> None:
        self._dict = None
        if _id is not None:
            samples = mongo_interface.get_samples(sample_ids=[ObjectId(_id)])
            if len(samples) == 1:
                self._dict = samples[0]
        elif name is not None:
            self._dict = {
                "name": name,
                "components": [],
                "properties": {
                    "paired_reads": {
                        "summary": {
                            "data": []
                        }
                    },
                    "sample_info": {
                        "summary": {
                            "emails": "",
                            "provided_species": "",
                            "comments": "",
                            "group": "",
                            "priority": ""
                        }
                    }
                },
                "report": {},
                "metadata": {}
            }

    def get(self, key: str) -> object:
        """Returns value of the associated key in object (dict)"""
        return self._dict[key]

    def set_name(self, name: str) -> None:
        """Sets name value in object (dict)"""
        self._dict["name"] = name

    # def set_components(self, components: list(Components)) -> None:
    #     self._dict["components"] = []
    #     if self.get("_id") != None:
    #         get_sample_component(sample_id=self._dict["_id"])
    #     for component in components:
    #         self._dict["components"].append({"_id": component.get_id(), "name": component.get_name(), "status": component.})

    def set_properties_paired_reads(self, paired_reads: Category) -> None:
        """Sets paired reads value in object (dict)"""
        if paired_reads.get_name() == "paired_reads":
            self._dict["properties"]["paired_reads"] = paired_reads.display()

    def set_properties_sample_info(self, sample_info: Category) -> None:
        """Sets sample_info value in object (dict)"""
        if sample_info.get_name() == "sample_info":
            self._dict["properties"]["sample_info"] = sample_info.display()

    def set_report(self, report: dict) -> None:
        """Sets report value in object (dict)"""
        self._dict["report"] = report

    def set_metadata(self, metadata: dict) -> None:
        """Sets metadata value in object (dict)"""
        self._dict["metadata"] = metadata

    def display(self) -> dict:
        """Return a copy of the object (dict)"""
        return self._dict.copy()

    def save(self) -> None:
        """Saves the in object (dict) to the DB"""
        self._dict = mongo_interface.dump_sample_info(self._dict)


class Run:
    """
    Run Object for saving in bifrost DB document "runs"
    A Run object is a collection of samples primarily for the purpose of quality control.
    Issues of samples that lack basic information are also saved. Issues include:
        - Duplicate sample IDs
        - Modified sample IDs
        - Unusued sample IDs
        - Sample IDs with no reads
    When Run.save() is ran the values are committed to the DB. 
    Run names are unique in the DB so an error is returned for duplicate values.
    """
    def __init__(self, _id: str = None, name: str = None) -> None:
        self._dict = None
        if _id is not None:
            runs = mongo_interface.get_runs(run_id=ObjectId(_id))
            if len(runs) == 1:
                self._dict = runs[0]
        elif name is not None:
            runs = mongo_interface.get_runs(names=[name])
            if len(runs) == 1:
                self._dict = runs[0]
            elif not runs:
                self._dict = {
                    "name": name,
                    "type": "default",
                    "path": None,
                    "samples": [],
                    "issues": {
                        "duplicate_samples": [],
                        "modified_samples": [],
                        "unused_files": [],
                        "samples_without_reads": [],
                        "samples_without_metadata": []
                    },
                    "Comments": ""
                }

    def get(self, key: str) -> object:
        """Returns value of the associated key in object (dict)"""
        return self._dict[key]

    def set_name(self, name: str) -> None:
        """Sets name value in object (dict)"""
        self._dict["name"] = name

    def set_type(self, type_: str) -> None:
        """Sets type value in object (dict)"""
        self._dict["type"] = type_

    def set_path(self, path: str) -> None:
        """Sets path value in object (dict)"""
        self._dict["path"] = path

    def set_samples(self, samples: List[Sample]) -> None:
        """Sets samples value in object (dict) by clearing out values then setting to new values"""
        self._dict["samples"] = []
        for sample in samples:
            self._dict["samples"].append({"_id": sample.get("_id"), "name": sample.get("name")})

    def set_issues(self, duplicate_samples: List[str], modified_samples: List[str], unused_files: List[str], samples_without_reads: List[str], samples_without_metadata: List[str]) -> None:
        """Sets issue(s) value in object (dict), requires all issues to be passed in"""
        self._dict["issues"]["duplicate_samples"] = duplicate_samples
        self._dict["issues"]["modified_samples"] = modified_samples
        self._dict["issues"]["unused_files"] = unused_files
        self._dict["issues"]["samples_without_reads"] = samples_without_reads
        self._dict["issues"]["samples_without_metadata"] = samples_without_metadata

    def set_comments(self, comments: str) -> None:
        """Sets comments value in object (dict)"""
        self._dict["Comments"] = comments

    def display(self) -> dict:
        """Return a copy of the object (dict)"""
        return self._dict.copy()

    def save(self) -> None:
        """Saves the in object (dict) to the DB"""
        self._dict = mongo_interface.dump_run_info(self._dict)


# class SampleComponent:
#     def __init__(self, sample_id: str, component_id: str) -> None:
#         component = Component().load(component_id)
#         sample = Sample().load(sample_id)
#         categories = list[Category] #summary output for each category it belongs to, under properties and report, results would be cumulative,
#                                     # Category/Property object would have a summary and report, component id is uneccesary as it's from the same component
#         _dict = {
#             "_id": None,
#             "component": {
#                 "name": 
#             }
#         }

class SampleComponentObj:
    def __init__(self, sample_id, component_id, path=None):
        self.sample_id = sample_id
        self.component_id = component_id
        self.path = path
        self.sample_db = get_sample(sample_id=self.sample_id)
        self.component_db = get_component(component_id=self.component_id)
        self.sample_component_db = get_sample_component(sample_id=self.sample_id, component_id=self.component_id)
        if self.sample_component_db == None:
            save_sample_component({
                "sample": {"_id": self.sample_db["_id"], "name": self.sample_db["name"]},
                "component": {"_id": self.component_db["_id"], "name": self.component_db["display_name"]},
                "path": path,
                "status": "initialized"
            })
            self.sample_component_db = get_sample_component(sample_id=self.sample_id, component_id=self.component_id)
        self.initialized()

    def load(self):
        #HACK: NOTE: Changed name to display name, obviously more needs to change. Leaving this here as a temp fix as the main fix is refactoring this object
        return (self.sample_db["name"], self.component_db["display_name"], self.component_db["install"], self.component_db["options"], self.component_db["resources"])

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
        paired_reads = self.get_sample_properties_by_category("paired_reads")
        if "data" in paired_reads:
            return (paired_reads["data"][0], paired_reads["data"][1])
        else:
            return ("/dev/null", "/dev/null")

    def get_resources(self):
        return self.component_db["resources"]

    def get_options(self):
        return self.component_db["options"]

    def check_requirements(self, output_file="requirements_met", log=None):
        no_failures = True
        if self.component_db["requirements"] is not None:
            requirements = pandas.json_normalize(self.component_db["requirements"], sep=".").to_dict(orient='records')[0]  # a little loaded of a line, get requirements from component_db, use the pandas json function to turn it into a 2d dataframe, then convert that to a dict of known depth 2, 0 is for our 1 and only sheet
            for requirement in requirements:
                category = requirement.split(".")[0]
                if category == "sample":
                    field = requirement.split(".")[1:]
                    expected_value = requirements[requirement]
                    if not self.requirement_met(self.sample_db, field, expected_value, log):
                        no_failures = False
                elif category == "component":
                    component_to_check = requirement.split(".")[1]
                    field = requirement.split(".")[2:]
                    expected_value = requirements[requirement]
                    s_c_db = get_sample_component(sample_id=self.sample_id,
                                              component_name=component_to_check)
                    if not self.requirement_met(s_c_db, field, expected_value, log):
                        no_failures = False
                else:
                    no_failures = False
        if no_failures:
            open(os.path.join(self.component_db["display_name"], output_file), "w+").close()
        else:
            self.requirements_not_met()

    def get_current_status(self):
        return self.sample_component_db["status"]

    def update_status_in_sample_and_sample_component(self, status):
        self.sample_component_db["status"] = status 
        status_set = False
        # TODO: this code should be refactored out (next 2 lines) as it's because a sample isn't initiated which should be solved by a sampleObj
        if "components" not in self.sample_db:
            self.sample_db["components"] = []
        for component in self.sample_db["components"]:
            if component["_id"] == self.component_db["_id"]:
                component["status"] = status
                status_set = True
        if not status_set:
            self.sample_db["components"].append({"_id": self.component_db["_id"], "name":self.component_db["display_name"], "status":status})
        self.save()

    def requirements_not_met(self):
        sys.stdout.write("Workflow stopped due to requirements\n")
        self.update_status_in_sample_and_sample_component("Requirements not met")

    def initialized(self):
        sys.stdout.write("Workflow initialized\n")
        self.update_status_in_sample_and_sample_component("Initialized")

    def queued(self):
        sys.stdout.write("Workflow queue'd\n")
        self.update_status_in_sample_and_sample_component("Queued")

    def started(self):
        sys.stdout.write("Workflow processing\n")
        self.update_status_in_sample_and_sample_component("Running")

    def success(self):
        sys.stdout.write("Workflow complete\n")
        if self.get_current_status() == "Running":
            self.update_status_in_sample_and_sample_component("Success")

    def failure(self):
        sys.stdout.write("Workflow error\n")
        if self.get_current_status() == "Running":
            self.update_status_in_sample_and_sample_component("Failure")

    def start_rule(self, rule_name, log=None):
        self.write_log_out(log, "{} has started\n".format(rule_name))
        return (self.component_db["display_name"], self.component_db["options"], self.component_db["resources"])

    def rule_run_cmd(self, command, log):
        self.write_log_out(log, "Running:{}\n".format(command))
        command_log_out, command_log_err = subprocess.Popen(command, shell=True).communicate()
        self.write_log_out(log, command_log_out)
        self.write_log_err(log, command_log_err)

    def end_rule(self, rule_name, log=None):
        self.write_log_out(log, "{} has finished\n".format(rule_name))

    def start_data_dump(self, log=None):
        if "properties" in self.sample_component_db:
            self.sample_component_db["properties"] = {
                "summary": {},
                "component": {
                    "_id": self.component_db["_id"]
                }
            }
        else:
            self.sample_component_db["properties"] = {
                "summary": {},
                "component": {
                    "_id": self.component_db["_id"]
                }
            }
        self.sample_component_db["results"] = {}
        if self.component_db["db_values_changes"]["sample"].get("report", None) is not None:
            # HACK: Right now we only support 1 category, thus the index 0 ref, when this is refactored we can support multiple. Config is already being adjusted to handle multiple
            self.sample_component_db["report"] = self.component_db["db_values_changes"]["sample"]["report"][self.component_db["category"][0]] 
        else:
            self.sample_component_db["report"] = {}
        self.write_log_out(log, "Starting datadump\n")
        self.save_files_to_sample_component(log)

    def save_files_to_sample_component(self, log=None):
        save_files_to_db(self.component_db["db_values_changes"]["files"], sample_component_id=self.sample_component_db["_id"])
        self.write_log_out(log, "Files saved: {}\n".format(",".join(self.component_db["db_values_changes"]["files"])))

    def get_component_name(self):
        return self.component_db["display_name"]

    def run_data_dump_on_function(self, data_extraction_function, log=None):
        (self.sample_component_db["properties"]["summary"], self.sample_component_db["results"]) = data_extraction_function(self)

    def end_data_dump(self, output_file="datadump_complete", generate_report_function=lambda x: None, log=None):
        # HACK: Right now we only support 1 category, thus the index 0 ref, when this is refactored we can support multiple. Config is already being adjusted to handle multiple
        self.sample_db["properties"][self.component_db["category"][0]] = self.sample_component_db["properties"]
        report_data = generate_report_function(self)
        if report_data is not None:
            self.sample_db["report"][self.component_db["category"][0]] = self.sample_component_db["report"]
            self.sample_db["report"][self.component_db["category"][0]]["data"] = report_data
            assert(type(self.sample_db["report"][self.component_db["category"][0]]["data"])==list)
        self.write_log_err(log, str(traceback.format_exc()))
        self.save()
        self.write_log_out(log, "sample {} saved\nsample_component {} saved\n".format(self.sample_db["_id"], self.sample_component_db["_id"]))
        open(os.path.join(self.component_db["display_name"], output_file), "w+").close()
        self.write_log_out(log, "Done datadump\n")
        self.success()

    def save(self):
        save_sample(self.sample_db)
        save_sample_component(self.sample_component_db)

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
        if content is not None:
            if log is not None:
                self.write_log(log.out_file, content)
            else:
                sys.stdout.write(content)

    def write_log_err(self, log, content):
        if content is not None:
            if log is not None:
                self.write_log(log.err_file, content)
            else:
                sys.stderr.write(content)

    def write_log(self, log_file, content):
        if content is not None:
            with open(log_file, "a+") as file_handle:
                file_handle.write(content)


def check_db_connection_exists():
    try:
        connection = mongo_interface.get_connection()
        connection.server_info()
        return True
    except:
        print(str(traceback.format_exc()))
        return False


def get_connection_info():
    connection = mongo_interface.get_connection()
    message = (
        f"Connected to:\n"
        f"    Database: {connection.get_database().name}\n"
        f"    Host: {':'.join([str(i) for i in connection.address])}\n"
    )
    return message


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


def get_sample_component(sample_component_id=None, sample_id=None, component_id=None, component_name=None):
    if sample_component_id is not None:
        sample_component_id = [ObjectId(sample_component_id)]
    if sample_id is not None:
        sample_id = [ObjectId(sample_id)]
    if component_id is not None:
        component_id = [ObjectId(component_id)]
    if component_name is not None:
        component_name = [component_name]
    return next(iter(mongo_interface.get_sample_components(sample_component_ids=sample_component_id, sample_ids=sample_id, component_ids=component_id, component_names=component_name)), None)


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
        load_file_from_db(f._id, save_to_path, subpath=True)

def recreate_yaml(collection, oid, path=None):
    if collection == "sample_components":
        obj = get_sample_component(oid)
        filename = "{}__{}.yaml".format(obj["sample"]["name"], obj["component"]["name"])
    elif collection == "samples":
        obj = get_sample(oid)
        filename = "sample.yaml"
    elif collection == "runs":
        obj = get_runs(oid)[0]
        filename = "run.yaml"
    else:
        raise ValueError("invalid collection")
    if path is None:
        path = "."
    path = os.path.join(path, filename)
    with open(path, "w") as file_handle:
        yaml.dump(obj, file_handle)


def recreate_s_c(sample_component_id, save_to_path):
    recreate_s_c_files(sample_component_id, save_to_path)
    recreate_yaml("sample_components", sample_component_id, save_to_path)

def recreate_sample(sample_id, save_to_path):
    sample = get_sample(sample_id)
    if sample is None:
        raise ValueError("Sample not found")
    name = sample["name"]
    path = os.path.join(save_to_path, name)
    if os.path.exists(path):
        raise FileExistsError("Sample directory already exists")
    os.makedirs(path)
    recreate_yaml("samples", sample_id, path)
    
    sample_components = get_sample_components(sample_ids=[sample_id])
    for sample_component in sample_components:
        recreate_s_c(str(sample_component["_id"]), path)

def recreate_run(run_id, save_to_path):
    runs = get_runs(run_id)
    if len(runs) == 0:
        raise ValueError("Run not found")
    run = runs[0]
    name = run["name"]
    path = os.path.join(save_to_path, name)
    if os.path.exists(path):
        raise FileExistsError("Run directory already exists")
    run_yaml_dir = os.path.join(path, "bifrost")
    os.makedirs(run_yaml_dir)
    recreate_yaml("runs", run_id, run_yaml_dir)
    for sample in run["samples"]:
        recreate_sample(str(sample["_id"]), path)
    #component.yamls?
    #create samples folder

def test():
    print("Hello")
    return (u'This is the bifrost lib')
