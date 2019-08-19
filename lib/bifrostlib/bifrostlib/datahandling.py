import os
import ruamel.yaml
from bson.objectid import ObjectId
from bson.int64 import Int64
from bifrostlib import mongo_interface
import pymongo
import traceback

# -----Deciding whether to keep this --------
# class DatadumpSampleComponentObj:
#     def __init__(self, folder, sample_file, component_file, sample_component_file, log):
#         self.folder = folder
#         self.sample_file = sample_file
#         self.component_file = component_file
#         self.sample_component_file = sample_component_file
#         self.log = log
#         self.load()
#
#     def load(self):
#         self.db_sample = load_sample(self.sample_file)
#         self.db_component = load_component(self.component_file)
#         self.db_sample_component = load_sample_component(self.sample_component_file)
#         self.db_sample_component["summary"] = self.db_sample_component.get("summary", {})
#         self.db_sample_component["results"] = self.db_sample_component.get("results", {})
#
#     def save(self):
#         save_sample(self.db_sample, self.sample_file)
#         save_sample_component(self.db_sample_component, self.sample_component_file)
#
#     def log_out(self, str(content)):
#         log(self.log.out_file, content)
#
#     def log_err(self, str(content)):
#         log(self.log.err_file, content)

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
    if not os.path.isfile(file_yaml):
        return {}
    with open(file_yaml, "r") as file_handle:
        return yaml.load(file_handle)


def save_sample(sample_dict, file_yaml):
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


def save_sample_component(sample_component_dict, file_yaml):
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


def sample_component_success(file_yaml, component):
    sample_component_dict = load_sample_component(file_yaml)
    if sample_component_dict["status"] != "Success":
        return False
    else:
        return True


def update_sample_component_success(file_yaml, component):
    sample_component_dict = load_sample_component(file_yaml)
    sample_component_dict["status"] = "Success"
    save_sample_component(sample_component_dict, file_yaml)


def update_sample_component_failure(file_yaml, component):
    sample_component_dict = load_sample_component(file_yaml)
    if sample_component_dict["status"] != "Requirements not met":
        sample_component_dict["status"] = "Failure"
    save_sample_component(sample_component_dict, file_yaml)


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


def datadump_template(extraction_callback, db, temp_data={}, key=None, file_path=None):
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


def get_samples(sample_ids=None, run_names=None,
                component_ids=None):
    if sample_ids is not None:
        sample_ids = [ObjectId(id) for id in sample_ids]
    if component_ids is not None:
        component_ids = [ObjectId(id) for id in component_ids]
    return mongo_interface.get_samples(sample_ids=sample_ids,
                                       run_names=run_names,
                                       component_ids=component_ids)


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


def get_sample_components(sample_component_ids=None,
                          sample_ids=None, component_names=None):
    # Should be smarter
    if sample_component_ids is not None:
        sample_component_ids = [ObjectId(id) for id in sample_component_ids]
    if sample_ids is not None:
        sample_ids = [ObjectId(id) for id in sample_ids]
    return mongo_interface.get_sample_components(
        sample_component_ids=sample_component_ids,
        sample_ids=sample_ids,
        component_names=component_names)


def post_sample_component(sample_component):
    return mongo_interface.dump_sample_component_info(sample_component)


def delete_sample_component(s_c_id=None, sample_id=None):
    if s_c_id is not None:
        s_c_id = ObjectId(s_c_id)
    if sample_id is not None:
        sample_id = ObjectId(sample_id)
    return mongo_interface.delete_sample_component(s_c_id, sample_id)


# /component

def get_components(component_ids=None):
    """
    Get components from db
    """
    if component_ids is not None:
        component_ids = list(map(ObjectId, component_ids))
    return mongo_interface.get_components(component_ids=component_ids)


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
