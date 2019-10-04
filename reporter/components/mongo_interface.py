import os
from datetime import datetime
import math
import re
import pymongo
import keys  # .gitgnored file
from bson.objectid import ObjectId
from bson.son import SON
import atexit

PAGESIZE = 10


def date_now():
    """
    Needed to keep the same date in python and mongo, as mongo rounds to millisecond
    """
    d = datetime.utcnow()
    return d.replace(microsecond=math.floor(d.microsecond/1000)*1000)

CONNECTION = None


def close_connection():
    global CONNECTION
    if CONNECTION is not None:
        CONNECTION.close()


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



def check_run_name(name):
    connection = get_connection()
    db = connection.get_database()
    # Fastest.
    run = db.runs.find({"name": name}).limit(1).count(True)
    return run is not 0


def get_run_list():
    connection = get_connection()
    db = connection.get_database()
    # Fastest.
    runs = list(db.runs.find( {},#{"type": "routine"}, #Leave in routine
                                {"name": 1,
                                "_id": 0,
                                "samples": 1}).sort([['metadata.created_at', pymongo.DESCENDING]]))
    return runs


def get_group_list(run_name=None):
    connection = get_connection()
    db = connection.get_database()
    if run_name is not None:
        run = db.runs.find_one(
            {"name": run_name},
            {
                "_id": 0,
                "samples._id": 1
            }
        )
        if run is None:
            run_samples = []
        else:
            run_samples = run["samples"]
        sample_ids = [s["_id"] for s in run_samples]
        groups = list(db.samples.aggregate([
            {
                "$match": {
                    "_id": {"$in": sample_ids},
                }
            },
            {
                "$group": {
                    "_id": "$properties.sample_info.summary.group",
                    "count": {"$sum": 1}
                }
            }
        ]))
    else:
        groups = list(db.samples.aggregate([
            {
                "$group": {
                    "_id": "$properties.sample_info.summary.group",
                    "count": { "$sum":1 }
                }
            }
        ]))

    return groups


def get_species_list(species_source, run_name=None):
    connection = get_connection()
    db = connection.get_database()
    if species_source == "provided":
        spe_field = "properties.sample_info.summary.provided_species"
    else:
        spe_field = "properties.species_detection.summary.detected_species"
    if run_name is not None:
        run = db.runs.find_one(
            {"name": run_name},
            {
                "_id": 0,
                "samples._id": 1
            }
        )
        if run is None:
            run_samples = []
        else:
            run_samples = run["samples"]
        sample_ids = [s["_id"] for s in run_samples]
        species = list(db.samples.aggregate([
            {
                "$match": {
                    "_id": {"$in": sample_ids}
                }
            },
            {
                "$group": {
                    "_id": "$" + spe_field,
                    "count": {"$sum": 1}
                }
            },
            {
                "$sort": {"_id": 1}
            }
        ]))
    else:
        species = list(db.samples.aggregate([
            {
                "$group": {
                    "_id": "$" + spe_field,
                    "count": {"$sum": 1}
                }
            },
            {
                "$sort": {"_id": 1}
            }
        ]))
    return species


def filter_qc(qc_list):
    if qc_list is None or len(qc_list) == 0:
        return None
    qc_query = []
    for elem in qc_list:
        if elem == "Not checked":
            qc_query.append({"$and": [
                {"properties.datafiles.summary.paired_reads": {"$exists": True}},
                {"properties.stamper.summary.stamp.value": {"$exists": False}}
            ]})
        elif elem == "fail:core facility":
            qc_query.append({"$or": [
                        {"properties.datafiles.summary.paired_reads": {"$exists": False}},
                        {"properties.stamper.summary.stamp.value": "fail:core facility"}
                    ]
                })
        else:
            qc_query.append({"properties.stamper.summary.stamp.value": elem})
    return {"$match": {"$and": qc_query}}


def filter(run_names=None,
           species=None, species_source="species", group=None,
           qc_list=None, samples=None, pagination=None,
           sample_names=None,
           projection=None):
    if species_source == "provided":
        spe_field = "properties.sample_info.summary.provided_species"
    elif species_source == "detected":
        spe_field = "properties.species_detection.summary.detected_species"
    else:
        spe_field = "properties.species_detection.summary.species"
    connection = get_connection()
    db = connection.get_database()
    query = []
    sample_set = set()
    if sample_names is not None and len(sample_names) != 0:
        sample_names_query = []
        for s_n in sample_names:
            if s_n.startswith("/") and s_n.endswith("/"):
                sample_names_query.append(re.compile(s_n[1:-1]))
            else:
                sample_names_query.append(s_n)
        query.append({"name": {"$in": sample_names_query}})
    if samples is not None and len(samples) != 0:
        sample_set = {ObjectId(id) for id in samples}
        query.append({"_id": {"$in": list(sample_set)}})
    if run_names is not None and len(run_names) != 0:
        runs = list(db.runs.find(
            {"name": {"$in": run_names}},
            {
                "_id": 0,
                "samples._id": 1
            }
        ))
        if runs is None:
            run_sample_set = set()
        else:
            run_sample_set = {s["_id"] for run in runs for s in run['samples']}
    
        if len(sample_set):
            inter = run_sample_set.intersect(sample_set)
            query.append({"_id": {"$in": list(inter)}})
        else:
            query.append({"_id": {"$in": list(run_sample_set)}})
    if species is not None and len(species) != 0:
        

        if "Not classified" in species:
            query.append({"$or":
                [
                    {spe_field: None},
                    {spe_field: {"$in": species}},
                    {spe_field: {"$exists": False}}
                ]
            })
        else:
            query.append({spe_field: {"$in": species}})
    if group is not None and len(group) != 0:
        if "Not defined" in group:
            query.append({"$or":
                [
                    {"properties.sample_info.summary.group": None},
                    {"properties.sample_info.summary.group": {"$in": group}},
                    {"properties.sample_info.summary.group": {"$exists": False}}
                ]
            })
        else:
            query.append(
                {"properties.sample_info.summary.group": {"$in": group}})

    if pagination is not None:
        p_limit = pagination['page_size']
        p_skip = pagination['page_size'] * pagination['current_page']
    else:
        p_limit = 1000
        p_skip = 0

    skip_limit_steps = [
        {"$skip": p_skip}, {"$limit": p_limit}
    ]

    qc_query = filter_qc(qc_list)

    if len(query) == 0:
        if qc_query is None:
            match_query = {}
        else:
            match_query = qc_query["$match"]
    else:
        if qc_query is None:
            match_query = {"$and": query}
        else:
            match_query = {"$and": query + qc_query["$match"]["$and"]}

    query_result = list(db.samples.find(
        match_query, projection).sort([('name', pymongo.ASCENDING)]).skip(p_skip).limit(p_limit))

    return query_result

def get_sample_runs(sample_ids):
    connection = get_connection()
    db = connection.get_database()
    return list(db.runs.find({"samples": {"$elemMatch": {"_id": {"$in": sample_ids}}}}))


def get_read_paths(sample_ids):
    connection = get_connection()
    db = connection.get_database()
    return list(db.samples.find({"_id": {"$in": list(map(lambda x: ObjectId(x), sample_ids))}}, {"reads": 1, "name": 1}))


def get_assemblies_paths(sample_ids):
    connection = get_connection()
    db = connection.get_database()
    return list(db.sample_components.find({
        "sample._id": {"$in": list(map(lambda x: ObjectId(x), sample_ids))},
        "component.name": "assemblatron"
    }, {"path": 1, "sample": 1}))



def get_species_QC_values(ncbi_species):
    connection = get_connection()
    db = connection.get_database('bifrost_species')
    species = db.species.find_one({"ncbi_species": ncbi_species}, {
                        "min_length": 1, "max_length": 1})
    if species is not None:
        return species
    species = db.species.find_one({"organism": ncbi_species}, {
        "min_length": 1, "max_length": 1})
    if species is not None:
        return species
    species = db.species.find_one({"group": ncbi_species}, {
        "min_length": 1, "max_length": 1})
    if species is not None:
        return species
    species = db.species.find_one({"organism": "default"}, {
        "min_length": 1, "max_length": 1})
    return species


def get_sample_QC_status(last_runs):
    connection = get_connection()
    db = connection.get_database()
    samples = [sample
                for run in last_runs
                for sample in run["samples"]]
        
    samples_full = db.samples.find({"_id": {"$in": list(map(lambda x: x["_id"], samples))}},
                                   {"properties.stamper": 1,
                                    "properties.datafiles": 1,
                                    "name": 1})
    samples_by_ids = {str(s["_id"]): s for s in samples_full}

    samples_runs_qc = {}
    for sample in samples:
        sample_dict = {}
        if str(sample["_id"]) not in samples_by_ids:
            print("Missing sample from DB: " + str(sample["_id"]))
            continue
        name = samples_by_ids[str(sample["_id"])]["name"]
        for run in last_runs:
            for run_sample in run["samples"]:
                if name == samples_by_ids[str(run_sample["_id"])]["name"]:
                    sample_db = samples_by_ids.get(str(run_sample["_id"]), None)
                    if sample_db is not None:
                        qc_val = sample_db.get("properties", {}).get("stamper", {}).get(
                            "summary", {}).get("stamp", {}).get("value", "N/A")
                        reads = sample_db.get("properties", {}).get("datafiles", {}).get(
                            "summary", {}).get("paired_reads", [])

                        if qc_val == "N/A" and not reads:
                            qc_val = "CF(LF)"
                        expert_check = False
                        # if "supplying_lab_check" in stamps and "value" in stamps["supplying_lab_check"]:
                        #     qc_val = stamps["supplying_lab_check"]["value"]
                        #     expert_check = True

                        if qc_val == "supplying lab":
                            qc_val = "SL"
                        elif (qc_val == "core facility" or
                                qc_val == "resequence"):
                            qc_val = "CF"
                        elif qc_val == "OK" or qc_val == "accepted":
                            qc_val = "OK"

                        if expert_check:
                            qc_val += "*"
                        sample_dict[run["name"]] = qc_val
        samples_runs_qc[name] = sample_dict
    return samples_runs_qc


def get_last_runs(run, n, runtype):
    connection = get_connection()
    db = connection.get_database()

    run = db.runs.find_one({"name": run})
    run_date = run.get("metadata", {}).get("created_at")

    if run_date is not None:
        if runtype is not None:
            query = {"metadata.created_at": {"$lte": run_date}, "type": runtype}
        else:
            query = {"metadata.created_at": {"$lte": run_date}}
    else:
        if runtype is not None:
            query = {"type": runtype}
        else:
            query = {}
    return list(db.runs.find(query, {"name": 1, "samples": 1}).sort([['metadata.created_at', pymongo.DESCENDING]]).limit(n))


def get_samples(sample_id_list, projection=None):
    connection = get_connection()
    db = connection.get_database()
    if projection is None:
        projection = {}
    return list(db.samples.find({"_id": {"$in": sample_id_list}}, projection))

def get_sample(sample_id):
    connection = get_connection()
    db = connection.get_database()
    return db.samples.find_one({"_id": sample_id})


def save_sample(data_dict):
    """COPIED FROM BIFROSTLIB Insert sample dict into mongodb.
    Return the dict with an _id element"""
    connection = get_connection()
    db = connection.get_database()
    samples_db = db.samples  # Collection name is samples
    if "_id" in data_dict:
        data_dict = samples_db.find_one_and_update(
            filter={"_id": data_dict["_id"]},
            update={"$set": data_dict},
            # return new doc if one is upserted
            return_document=pymongo.ReturnDocument.AFTER,
            upsert=False  # insert the document if it does not exist
        )
    else:
        data_dict = samples_db.find_one_and_update(
            filter=data_dict,
            update={"$setOnInsert": data_dict},
            # return new doc if one is upserted
            return_document=pymongo.ReturnDocument.AFTER,
            upsert=True  # insert the document if it does not exist
        )
    return data_dict


def save_sample_component(data_dict):
    """COPIED FROM BIFROSTLIB. Insert sample dict into mongodb.
    Return the dict with an _id element"""
    connection = get_connection()
    db = connection.get_database()
    sample_components_db = db.sample_components
    now = date_now()
    data_dict["metadata"] = data_dict.get("metadata", {'created_at': now})
    data_dict["metadata"]["updated_at"] = now
    if "_id" in data_dict:
        data_dict = sample_components_db.find_one_and_update(
            filter={"_id": data_dict["_id"]},
            update={"$set": data_dict},
            # return new doc if one is upserted
            return_document=pymongo.ReturnDocument.AFTER,
            # This might change in the future. It doesnt make much sense with our current system.
            upsert=True
            # Import relies on this to be true.
            # insert the document if it does not exist
        )
    else:
        search_fields = {
            "sample._id": data_dict["sample"]["_id"],
            "component._id": data_dict["component"]["_id"],
        }
        data_dict = sample_components_db.find_one_and_update(
            filter=search_fields,
            update={
                "$set": data_dict
            },
            # return new doc if one is upserted
            return_document=pymongo.ReturnDocument.AFTER,
            upsert=True  # insert the document if it does not exist
        )
    return data_dict


def save_component(data_dict):
    """COPIED FROM BIFROSTLIB. Insert sample dict into mongodb.
    Return the dict with an _id element"""
    connection = get_connection()
    db = connection.get_database()
    now = date_now()
    data_dict["metadata"] = data_dict.get("metadata", {'created_at': now})
    data_dict["metadata"]["updated_at"] = now
    components_db = db.components  # Collection name is samples
    if "_id" in data_dict:
        data_dict = components_db.find_one_and_update(
            filter={"_id": data_dict["_id"]},
            update={"$set": data_dict},
            # return new doc if one is upserted
            return_document=pymongo.ReturnDocument.AFTER,
            upsert=True  # This might change in the future # insert the document if it does not exist
        )
    else:
        data_dict = components_db.find_one_and_update(
            filter=data_dict,
            update={"$setOnInsert": data_dict},
            # return new doc if one is upserted
            return_document=pymongo.ReturnDocument.AFTER,
            upsert=True  # insert the document if it does not exist
        )
    return data_dict

def get_run(run_name):
    connection = get_connection()
    db = connection.get_database()
    return db.runs.find_one({"name": run_name})


def get_component(name=None, version=None):
    connection = get_connection()
    db = connection.get_database()
    query = {}
    if name is not None:
        query["name"] = name

    if version is not None:
        query["version"] = version
    return db.components.find_one(
        query)

def get_comment(run_id):
    connection = get_connection()
    db = connection.get_database()
    return db.runs.find_one(
        {"_id": run_id}, {"Comments": 1})

def set_comment(run_id, comment):
    connection = get_connection()
    db = connection.get_database()
    ret = db.runs.find_one_and_update({"_id": run_id}, {"$set": {"Comments": comment}})
    if ret != None:
        return 1
    else:
        return 0
