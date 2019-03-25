import pkg_resources
import ruamel.yaml
from bifrostlib import mongo_interface
from bson.objectid import ObjectId
import os


ObjectId.yaml_tag = u'!bson.objectid.ObjectId'
ObjectId.to_yaml = classmethod(
    lambda cls, representer, node: representer.represent_scalar(
        cls.yaml_tag, u'{}'.format(node)))
ObjectId.from_yaml = classmethod(
    lambda cls, constructor, node: cls(node.value))

yaml = ruamel.yaml.YAML(typ="safe")
yaml.default_flow_style = False
yaml.register_class(ObjectId)


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


def datadump_template(data_dict, component_folder,
                      file_path, extraction_callback):
    file_path_key = file_path.replace(".", "_")
    if os.path.isfile(os.path.join(component_folder, file_path)):
        data_dict["results"][file_path_key] = {}
        try:
            data_dict = extraction_callback(os.path.join(component_folder,
                                                         file_path),
                                            file_path_key, data_dict)
        except Exception as e:
            print(file_path, e)
            data_dict["results"][file_path_key]["status"] = "datadumper error"
        return data_dict


# /runs

def get_runs(names=None,sample_id=None):
    return mongo_interface.get_runs(names=names,
                                    sample_id=ObjectId(sample_id))

def delete_run(name=None):
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


def get_sample_components(sample_ids=None, component_names=None):
    # Should be smarter
    if sample_ids is not None:
        sample_ids = [ObjectId(id) for id in sample_ids]
    return mongo_interface.get_sample_components(sample_ids,
                                                 component_names)


def save_sample_component_to_db(sample_component):
    return mongo_interface.dump_sample_component_info(sample_component)


def delete_sample_component(s_c_id=None, sample_id=None):
    if s_c_id is not None:
        s_c_id = ObjectId(s_c_id)
    if sample_id is not None:
        sample_id = ObjectId(sample_id)
    return mongo_interface.delete_sample_component(s_c_id, sample_id)


# /component


def save_component_to_db(component):
    return mongo_interface.dump_component_info(component)


def delete_component(component_id):
    return mongo_interface.delete_component(ObjectId(component_id))


# /species
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


def test():
    print("Hello")
    return (u'This is the bifrost lib')
