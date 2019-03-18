import pymongo
import re
import os
import ruamel.yaml
import traceback
yaml = ruamel.yaml.YAML(typ="safe")
yaml.default_flow_style = False


def get_connection():
    # config_file = "run_config.yaml"
    # if not os.path.isfile("run_config.yaml"):
    #     config_file = "../run_config.yaml"
    # with open(config_file, "r") as file_handle:
    #     config = yaml.load(file_handle)
    mongo_db_key_location = os.getenv("BIFROST_DB_KEY", None)
    with open(mongo_db_key_location, "r") as mongo_db_key_location_handle:
        mongodb_url = mongo_db_key_location_handle.readline().strip()
    "Return mongodb connection"
    return pymongo.MongoClient(mongodb_url)


def get_species_connection():
    mongo_db_key_location = os.getenv("BIFROST_SPECIES_DB_KEY", None)
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
        with get_species_connection() as connection:
            db = connection.get_database()  # Database name is ngs_runs
            species_db = db.species  # Collection name is samples
            result = species_db.find_one({"ncbi_species": ncbi_species_name}, {"mlst_species": 1, "_id": 0})
            if result is not None and "mlst_species" in result:
                return result["mlst_species"]
            else:
                return None
    except Exception as e:
        print(traceback.format_exc())
        return None


def query_ncbi_species(species_entry):
    if species_entry is None:
        return None
    try:
        species_entry = re.compile(species_entry)
        with get_species_connection() as connection:
            db = connection.get_database()  # Database name is ngs_runs
            species_db = db.species  # Collection name is samples
            result = species_db.find_one({"organism": species_entry}, {"ncbi_species": 1, "_id": 0})
            group_result = species_db.find_one({"group": species_entry}, {"ncbi_species": 1, "_id": 0})
            if result is not None:
                return result["ncbi_species"]
            elif group_result is not None:
                return group_result["ncbi_species"]
            else:
                return None
    except Exception as e:
        print(traceback.format_exc())
        return None


def query_species(ncbi_species_name):
    try:
        with get_species_connection() as connection:
            db = connection.get_database()  # Database name is ngs_runs
            species_db = db.species
            if ncbi_species_name is not None and species_db.find({'ncbi_species': ncbi_species_name}).count() > 0:
                return species_db.find_one({"ncbi_species": ncbi_species_name})
            else:
                return species_db.find_one({"organism": "default"})
    except Exception as e:
        print(traceback.format_exc())
        return None


def load_run(name=None):
    try:
        with get_connection() as connection:
            db = connection.get_database()
            if name is not None:
                return db.runs.find_one({"name": name, "type": "routine"})
            else:
                return db.runs.find_one(
                    {"$query": {"type": "routine"}, "$orderby": {"_id": -1}})
    except Exception as e:
        print(traceback.format_exc())
        return None


def load_sample(sample_id):
    try:
        with get_connection() as connection:
            db = connection.get_database()
            return db.samples.find_one({"_id": sample_id})
    except Exception as e:
        print(traceback.format_exc())
        return None


def load_samples(sample_ids):
    try:
        with get_connection() as connection:
            db = connection.get_database()
            return list(db.samples.find({"_id": {"$in": sample_ids}}))
    except Exception as e:
        print(traceback.format_exc())
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
        print(traceback.format_exc())
        return None


def load_samples_from_runs(run_ids=None, names=None):
    try:
        with get_connection() as connection:
            db = connection.get_database()
            if run_ids is not None:
                return list(db.runs.find({"_id": {"$in": run_ids}}, {"samples": 1}))
            elif names is not None:
                return list(db.runs.find({"name": {"$in": names}}, {"samples": 1}))
            else:
                return [db.runs.find_one({"$query": {"type": "routine"}, "$orderby": {"_id": -1}})]
    except Exception as e:
        print(traceback.format_exc())
        return None


def load_all_samples():
    try:
        with get_connection() as connection:
            db = connection.get_database()
            runs = list(db.runs.find({}, {"samples": 1}))
            sample_ids = set()
            for run in runs:
                sample_ids.update(list(map(lambda x: x["_id"],
                                           run["samples"])))

        return list(sample_ids)

    except Exception as e:
        print(traceback.format_exc())
        return None


def delete_sample_component(s_c_id=None, sample_id=None,
                            s_c_id_list=None, sample_id_list=None):
    """
    DELETE sample component from database. Parameters are AND connected.
    Returns number of deleted entries
    """
    query = []
    if s_c_id is not None:
        query.append({"_id": s_c_id})
    if sample_id is not None:
        query.append({"sample._id": sample_id})
    if s_c_id_list is not None:
        query.append({"_id": {"$in": s_c_id_list}})
    if sample_id_list is not None:
        query.append({"sample._id": {"$in": sample_id_list}})
    try:
        with get_connection() as connection:
            db = connection.get_database()
            result = db.sample_components.delete_many({"$and": query})
            return result.deleted_count
    except Exception as e:
        print(traceback.format_exc())
        return None


def delete_sample_from_runs(sample_id=None):
    """
    Delete sample from runs. Returns number of runs this sample was deleted
    from.
    """
    update_count = 0
    try:
        with get_connection() as connection:
            db = connection.get_database()

            runs = db.runs.find({"samples._id": sample_id})
            for run in runs:
                new_samples = [sample
                               for sample in run["samples"]
                               if sample["_id"] != sample_id]
                run["samples"] = new_samples
                result = db.runs.update_one({"_id": run["_id"]}, run)
                update_count += result.modified_count
        return update_count
    except Exception as e:
        print(traceback.format_exc())
        return None


def delete_sample(sample_id):
    try:
        with get_connection() as connection:
            db = connection.get_database()

            result = db.samples.delete_one({"_id": sample_id})
        return result.deleted_count
    except Exception as e:
        print(traceback.format_exc())
        return None
