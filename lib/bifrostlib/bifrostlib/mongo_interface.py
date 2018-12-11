import pymongo
import re
import os
import ruamel.yaml
yaml = ruamel.yaml.YAML(typ="safe")
yaml.default_flow_style = False


def get_connection():
    config_file = "run_config.yaml"
    if not os.path.isfile("run_config.yaml"):
        config_file = "../run_config.yaml"
    with open(config_file, "r") as file_handle:
        config = yaml.load(file_handle)
    mongo_db_key_location = os.path.join(os.path.dirname(__file__), config["mongo_db_key_location"])
    with open(mongo_db_key_location, "r") as mongo_db_key_location_handle:
        mongodb_url = mongo_db_key_location_handle.readline().strip()
    "Return mongodb connection"
    return pymongo.MongoClient(mongodb_url)


def dump_run_info(data_dict):
    """Insert sample dict into mongodb.
    Return the dict with an _id element"""
    with get_connection() as connection:
        db = connection.get_database()
        runs_db = db.runs  # Collection name is samples
        if "_id" in data_dict:
            data_dict = runs_db.find_one_and_update(
                filter={"_id": data_dict["_id"]},
                update={"$set": data_dict},
                return_document=pymongo.ReturnDocument.AFTER,  # return new doc if one is upserted
                upsert=False  # insert the document if it does not exist
            )
        else:
            data_dict = runs_db.find_one_and_update(
                filter=data_dict,
                update={"$setOnInsert": data_dict},
                return_document=pymongo.ReturnDocument.AFTER,  # return new doc if one is upserted
                upsert=True  # insert the document if it does not exist
            )
        return data_dict


def dump_sample_info(data_dict):
    """Insert sample dict into mongodb.
    Return the dict with an _id element"""
    with get_connection() as connection:
        db = connection.get_database()
        samples_db = db.samples  # Collection name is samples
        if "_id" in data_dict:
            data_dict = samples_db.find_one_and_update(
                filter={"_id": data_dict["_id"]},
                update={"$set": data_dict},
                return_document=pymongo.ReturnDocument.AFTER,  # return new doc if one is upserted
                upsert=False  # insert the document if it does not exist
            )
        else:
            data_dict = samples_db.find_one_and_update(
                filter=data_dict,
                update={"$setOnInsert": data_dict},
                return_document=pymongo.ReturnDocument.AFTER,  # return new doc if one is upserted
                upsert=True  # insert the document if it does not exist
            )
        return data_dict


def dump_component_info(data_dict):
    """Insert sample dict into mongodb.
    Return the dict with an _id element"""
    with get_connection() as connection:
        db = connection.get_database()
        components_db = db.components  # Collection name is samples
        if "_id" in data_dict:
            data_dict = components_db.find_one_and_update(
                filter={"_id": data_dict["_id"]},
                update={"$set": data_dict},
                return_document=pymongo.ReturnDocument.AFTER,  # return new doc if one is upserted
                upsert=False  # insert the document if it does not exist
            )
        else:
            data_dict = components_db.find_one_and_update(
                filter=data_dict,
                update={"$setOnInsert": data_dict},
                return_document=pymongo.ReturnDocument.AFTER,  # return new doc if one is upserted
                upsert=True  # insert the document if it does not exist
            )
        return data_dict


def dump_sample_component_info(data_dict):
    """Insert sample dict into mongodb.
    Return the dict with an _id element"""
    with get_connection() as connection:
        db = connection.get_database()
        sample_components_db = db.sample_components  # Collection name is samples
        if "_id" in data_dict:
            data_dict = sample_components_db.find_one_and_update(
                filter={"_id": data_dict["_id"]},
                update={"$set": data_dict},
                return_document=pymongo.ReturnDocument.AFTER,  # return new doc if one is upserted
                upsert=False  # insert the document if it does not exist
            )
        else:
            data_dict = sample_components_db.find_one_and_update(
                filter=data_dict,
                update={"$setOnInsert": data_dict},
                return_document=pymongo.ReturnDocument.AFTER,  # return new doc if one is upserted
                upsert=True  # insert the document if it does not exist
            )
        return data_dict


def query_mlst_species(ncbi_species_name):
    if ncbi_species_name is None:
        return None
    try:
        ncbi_species_name = re.compile(ncbi_species_name)
        with get_connection() as connection:
            db = connection.get_database()  # Database name is ngs_runs
            species_db = db.species  # Collection name is samples
            result = species_db.find_one({"ncbi_species": ncbi_species_name}, {"mlst_species": 1, "_id": 0})
            if result is not None:
                return result["mlst_species"]
            else:
                return None
    except Exception as e:
        return None


def query_ncbi_species(species_entry):
    if species_entry is None:
        return None
    try:
        species_entry = re.compile(species_entry)
        with get_connection() as connection:
            db = connection.get_database()  # Database name is ngs_runs
            species_db = db.species  # Collection name is samples
            result = species_db.find_one({"organism": species_entry}, {"ncbi_species": 1, "_id": 0})
            if result["ncbi_species"] is not None:
                return result["ncbi_species"]
            else:
                # Returns real group name, instead of species_entry in case it contains
                in_group = species_db.find_one({"group": species_entry}, {"group": 1, "_id": 0})
                if in_group["group"] is not None:
                    return in_group["group"]
                else:
                    return None
    except Exception as e:
        return None

def query_species(ncbi_species_name):
    try:
        with get_connection() as connection:
            db = connection.get_database()  # Database name is ngs_runs
            species_db = db.species
            if ncbi_species_name is not None and species_db.find({'ncbi_species': ncbi_species_name}).count() > 0:
                return species_db.find_one({"ncbi_species": ncbi_species_name})
            else:
                return species_db.find_one({"organism": "default"})
    except Exception as e:
        return None

def load_sample(sample_id):
    try:
        with get_connection() as connection:
            db = connection.get_database()
            return db.samples.find_one({"_id": sample_id})
    except Exception as e:
        return None

def load_last_sample_component(sample_id, component_name):
    """Loads most recent sample component for a sample"""
    try:
        with get_connection() as connection:
            db = connection.get_database()
            result = list(db.sample_components.find({"sample._id": sample_id, "component.name": component_name}).sort([("setup_date", -1)]).limit(1))
            if len(result) > 0:
                return result[0]
            else:
                return None

    except Exception as e:
        return None


def load_samples_from_runs(run_ids):
    try:
        with get_connection() as connection:
            db = connection.get_database()
            if type(run_ids) == list:
                return list(db.runs.find({"_id": {"$in": run_ids}}, {"samples": 1}))
            else:
                return list(db.runs.find_one({"_id": run_ids}, {"samples": 1}))

    except Exception as e:
        return None

def load_all_samples():
    try:
        with get_connection() as connection:
            db = connection.get_database()
            runs = list(db.runs.find({}, {"samples": 1}))
            sample_ids = set()
            for run in runs:
                sample_ids.update(list(map(lambda x:x["_id"], run["samples"])))

        return list(sample_ids)

    except Exception as e:
        return None
