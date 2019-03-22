import pymongo
import re
import os
import ruamel.yaml
import traceback
import atexit
yaml = ruamel.yaml.YAML(typ="safe")
yaml.default_flow_style = False

CONNECTION = None
SPECIES_CONNECTION = None


def close_connection():
    global CONNECTION
    global SPECIES_CONNECTION
    if CONNECTION is not None:
        CONNECTION.close()
    if SPECIES_CONNECTION is not None:
        SPECIES_CONNECTION.close()

atexit.register(close_connection)


def get_connection():
    global CONNECTION
    if CONNECTION is not None:
        return CONNECTION
    else:
        mongo_db_key_location = os.getenv("BIFROST_DB_KEY", None)
        with open(mongo_db_key_location, "r") as mongo_db_key_location_handle:
            mongodb_url = mongo_db_key_location_handle.readline().strip()
        "Return mongodb connection"
        CONNECTION = pymongo.MongoClient(mongodb_url)
        return CONNECTION


def get_species_connection():
    global SPECIES_CONNECTION
    if SPECIES_CONNECTION is not None:
        return SPECIES_CONNECTION
    else:
        mongo_db_key_location = os.getenv("BIFROST_SPECIES_DB_KEY", None)
        with open(mongo_db_key_location, "r") as mongo_db_key_location_handle:
            mongodb_url = mongo_db_key_location_handle.readline().strip()
        "Return mongodb connection"
        SPECIES_CONNECTION = pymongo.MongoClient(mongodb_url)
        return SPECIES_CONNECTION


def dump_run_info(data_dict):
    """Insert sample dict into mongodb.
    Return the dict with an _id element"""
    connection = get_connection()
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
    connection = get_connection()
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
    connection = get_connection()
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
    connection = get_connection()
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


# Should call get_species
def query_mlst_species(ncbi_species_name):
    if ncbi_species_name is None:
        return None
    try:
        ncbi_species_name = re.compile(ncbi_species_name)
        connection = get_species_connection()
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


# Should call get_species
def query_ncbi_species(species_entry):
    if species_entry is None:
        return None
    try:
        species_entry = re.compile(species_entry)
        connection = get_species_connection()
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


# Should be renamed to get_species
def query_species(ncbi_species_name):
    try:
        connection = get_species_connection()
        db = connection.get_database()  # Database name is ngs_runs
        species_db = db.species
        if ncbi_species_name is not None and species_db.find({'ncbi_species': ncbi_species_name}).count() > 0:
            return species_db.find_one({"ncbi_species": ncbi_species_name})
        else:
            return species_db.find_one({"organism": "default"})
    except Exception as e:
        print(traceback.format_exc())
        return None


# DEPRECATED
def load_run(**kwargs):
    get_runs(**kwargs)


def get_runs(names=None, size=0):

    size = min(1000, size)
    size = max(-1000, size)

    query = {"type": "routine"}

    try:
        connection = get_connection()
        db = connection.get_database()
        if name is not None:
            query["name"] = name
        return list(db.runs.find(
            {"$query": query, "$orderby": {"_id": -1}}).limit(size))
    except Exception as e:
        print(traceback.format_exc())
        return None


def get_samples(sample_ids=None, run_names=None):
    # Uses AND operand

    query = []
    if sample_ids is not None:
        query.append({"_id": {"$in": sample_ids}})
    if run_names is not None:
        run = db.runs.find_one({"name": {"$in": run_names}}, {"samples._id": 1})
        if run is not None:
            run_sample_ids = [s["_id"] for s in run["samples"]]
        query.append({"_id": {"$in": run_sample_ids}})

    try:
        connection = get_connection()
        db = connection.get_database()
        return list(db.samples.find({"$and": query}))
    except Exception as e:
        print(traceback.format_exc())
        return None


def get_sample_components(sample_ids=None,
                          component_names=None,
                          size=1):
    """Loads most recent sample component for a sample"""
    size = min(1000, size)
    size = max(-1000, size)

    query = []
    if sample_ids is not None:
        query.append({"sample._id": {"$in": sample_ids}})
    if component_names is not None:
        query.append({"component.name": {"$in": component_names}})
    try:
        print(query)
        connection = get_connection()
        db = connection.get_database()
        result = list(db.sample_components.find({"$and": query})
                                          .sort([("setup_date", -1)])
                                          .limit(size))

    except Exception as e:
        print(traceback.format_exc())
        return None


# Should call get_runs
def load_samples_from_runs(run_ids=None, names=None):
    try:
        connection = get_connection()
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


# Should call get_samples
def load_all_samples():
    try:
        connection = get_connection()
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
        connection = get_connection()
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
        connection = get_connection()
        db = connection.get_database()

        runs = db.runs.find({"samples._id": sample_id})
        for run in runs:
            new_samples = [sample
                            for sample in run["samples"]
                            if sample["_id"] != sample_id]
            run["samples"] = new_samples
            result = db.runs.replace_one({"_id": run["_id"]}, run)
            update_count += result.modified_count
        return update_count
    except Exception as e:
        print(traceback.format_exc())
        return None


def delete_sample(sample_id):
    try:
        connection = get_connection()
        db = connection.get_database()

        result = db.samples.delete_one({"_id": sample_id})
        return result.deleted_count
    except Exception as e:
        print(traceback.format_exc())
        return None
