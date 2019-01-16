import pkg_resources
import ruamel.yaml
from bifrostlib import mongo_interface
import bson
import os


bson.objectid.ObjectId.yaml_tag = u'!bson.objectid.ObjectId'
bson.objectid.ObjectId.to_yaml = classmethod(lambda cls, representer, node: representer.represent_scalar(cls.yaml_tag, u'{}'.format(node)))
bson.objectid.ObjectId.from_yaml = classmethod(lambda cls, constructor, node: cls(node.value))

yaml = ruamel.yaml.YAML(typ="safe")
yaml.default_flow_style = False
yaml.register_class(bson.objectid.ObjectId)


def connect_or_initialize_and_connect(config):
    mongo_interface.set_connection_key_location(config)


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
    if len(file_yaml.split("/")) >= 2 and os.path.isdir("/".join(file_yaml.split("/")[:-1])) is False:
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
    if len(file_yaml.split("/")) >= 2 and os.path.isdir("/".join(file_yaml.split("/")[:-1])) is False:
        os.mkdir("/".join(file_yaml.split("/")[:-1]))
    if config["use_mongodb"]:
        sample_component_dict = mongo_interface.dump_sample_component_info(sample_component_dict)
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
    sample_component_dict["path"] = os.path.join(os.getcwd(), component)
    save_sample_component(sample_component_dict, file_yaml)


def update_sample_component_failure(file_yaml, component):
    sample_component_dict = load_sample_component(file_yaml)
    if sample_component_dict["status"] != "Requirements not met":
        sample_component_dict["status"] = "Failure"
    sample_component_dict["path"] = os.path.join(os.getcwd(), component)
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


def get_mlst_species_DB(file_yaml):
    with open(file_yaml, "r") as file_handle:
        sample = yaml.load(file_handle)
        return mongo_interface.query_mlst_species(sample["properties"].get("species", None))


def load_species(ncbi_species):
    return mongo_interface.query_species(ncbi_species)


def get_ncbi_species(species_entry):
    return mongo_interface.query_ncbi_species(species_entry)


def datadump_template(data_dict, component_folder, file_path, extraction_callback):
    file_path_key = file_path.replace(".", "_")
    if os.path.isfile(os.path.join(component_folder, file_path)):
        data_dict["results"][file_path_key] = {}
        try:
            data_dict = extraction_callback(os.path.join(component_folder, file_path), file_path_key, data_dict)
        except Exception as e:
            print(file_path, e)
            data_dict["results"][file_path_key]["status"] = "datadumper error"
        return data_dict

def load_run_from_db(name=None):
    return mongo_interface.load_run(name)

def load_sample_from_db(sample_id):
    return mongo_interface.load_sample(bson.objectid.ObjectId(sample_id))


def load_samples_from_db(sample_ids):
    sample_ids = [bson.objectid.ObjectId(id) for id in sample_ids]
    return mongo_interface.load_samples(sample_ids)

def save_sample_to_db(sample):
    return mongo_interface.dump_sample_info(sample)


def save_sample_component_to_db(sample_component):
    return mongo_interface.dump_sample_component_info(sample_component)


def save_component_to_db(component):
    return mongo_interface.dump_component_info(component)


def load_last_sample_component(sample_id, component_name):
    return mongo_interface.load_last_sample_component(bson.objectid.ObjectId(sample_id), component_name)


def load_samples_from_runs(run_ids=None, names=None):
    if run_ids is not None:
        if type(run_ids) == str:
            run_ids = [run_ids]
        run_ids = [bson.objectid.ObjectId(id) for id in run_ids]
    res = mongo_interface.load_samples_from_runs(run_ids, names)
    sample_ids = set()  # avoid dupes
    for run in res:
        sample_ids.update(list(map(lambda x: x["_id"], run["samples"])))
    return list(sample_ids)


def load_all_samples():
    return mongo_interface.load_all_samples()


def test():
    print("Hello")
    return (u'This is the bifrost lib')
