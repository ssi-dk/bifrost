import os
import pymongo
import keys  # .gitgnored file
from bson.objectid import ObjectId
from bson.son import SON
import atexit


PAGESIZE = 25

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
                                "samples": 1}).sort([['name', pymongo.DESCENDING]]))
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
                    "_id": "$sample_sheet.group",
                    "count": {"$sum": 1}
                }
            }
        ]))
    else:
        groups = list(db.samples.aggregate([
            {
                "$group": {
                    "_id": "$sample_sheet.group",
                    "count": { "$sum":1 }
                }
            }
        ]))

    return groups


def get_species_list(species_source, run_name=None):
    connection = get_connection()
    db = connection.get_database()
    if species_source == "provided":
        spe_field = "properties.provided_species"
    else:
        spe_field = "properties.detected_species"
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


def get_qc_list(run_name=None):
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
        qcs = list(db.sample_components.aggregate([
            {
                "$match": {
                    "sample._id": {"$in": sample_ids},
                    #"status": "Success",
                    "component.name": "ssi_stamper"
                }
            },
            {"$sort": {"sample._id": 1, "_id": 1}},
            {
                "$group": {
                    "_id": "$sample._id",
                    "action": {"$last": "$summary.assemblatron:action"}
                }
            },
            {
                "$group": {
                    "_id": "$action",
                    "count": {"$sum": 1}
                }
            },
            {
                "$sort": {"_id": 1}
            }
        ]))
    else:
        runs = list(db.runs.find({}, #{"type": "routine"},
                                 {"samples": 1}))
        sample_ids = set()
        for run in runs:
            for sample in run["samples"]:
                sample_ids.add(sample["_id"])
        sample_list = list(sample_ids)
        qcs = list(db.samples.aggregate([
            {
                "$match": {
                    "_id": {"$in": sample_list}
                }
            },
            {
                "$lookup": {
                    "from": "sample_components",
                    "let": {"sample_id": "$_id"},
                    "pipeline": [
                        {"$match": {
                            "component.name": "ssi_stamper",
                            "summary.assemblatron:action" : {"$exists" : True}
                            }},
                        { "$match": {
                                "$expr": {"$eq": ["$sample._id", "$$sample_id"]}
                            }
                        },
                        {"$project": {"summary.assemblatron:action" : 1}},
                        {"$sort": {"_id": -1}},
                        {"$limit": 1}
                    ],
                    "as": "sample_components"
                }
            },
            {
                "$unwind": "$sample_components"
            },
            {
                "$group": {
                    "_id": "$sample_components.summary.assemblatron:action",
                    "count": {"$sum": 1}
                }
            },
            {
                "$sort": {"_id": 1}
            }
        ]))

    return qcs

def filter_qc(db, qc_list, query):
    qc_query = [{"sample_components.summary.assemblatron:action": {"$in": qc_list}} ]
    if "Not checked" in qc_list:
        qc_query += [
            {"$and": [
                {"reads.R1": {"$exists": True}},
                {"$or": [
                    {"sample_components": {"$size": 0}}
                    # Should probably uncomment after Undetermined check is finished
                    # {"sample_components.status": {"$ne": "Success"}}
                ]}
            ]}
        ]
    if "skipped" in qc_list:
        qc_query += [
            {"sample_components.status": "initialized"}
        ]
        # Uncomment when we implement skipped for Undetermined
        # qc_query += [
        #     {"sample_components.status": "skipped"}
        # ]
    if "core facility" in qc_list:
        qc_query += [
            {"reads.R1": {"$exists": False}}
        ]

    if len(qc_query) > 1:
        qc_query = {"$or": qc_query}
    else:
        qc_query = qc_query[0]


    result = db.samples.aggregate([
        {
            "$match": query,
        },
        {
            "$lookup": {
                "from": "sample_components",
                "let": {"sample_id": "$_id"},
                "pipeline": [
                    {
                        "$match": {
                            "$expr": {
                                "$and": [
                                    {"$eq" : ["$component.name", "ssi_stamper"]},
                                    {"$eq": ["$sample._id", "$$sample_id"]}
                                ]
                            }
                        }
                    }
                ],
                "as": "sample_components"
            }
        },
        {
            "$match": qc_query
            
        }
    ])
    return list(result)


def filter(projection=None, run_name=None,
           species=None, species_source="species", group=None, qc_list=None, samples=None):
    if qc_list == ["OK", "core facility", "supplying lab", "skipped", "Not checked"]:
        qc_list = None
    if species_source == "provided":
        spe_field = "properties.provided_species"
    elif species_source == "detected":
        spe_field = "properties.detected_species"
    else:
        spe_field = "properties.species"
    connection = get_connection()
    db = connection.get_database()
    query = []
    sample_set = set()
    if samples is not None:
        sample_set = {ObjectId(id) for id in samples}
        query.append({"_id": {"$in": list(sample_set)}})
    if run_name is not None and run_name != "":
        run = db.runs.find_one(
            {"name": run_name},
            {
                "_id": 0,
                "samples._id": 1
            }
        )
        if run is None:
            run_sample_set = set()
        else:
            run_sample_set = {s["_id"] for s in run['samples']}
    
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
                    {"sample_sheet.group": None},
                    {"sample_sheet.group": {"$in": group}},
                    {"sample_sheet.group": {"$exists": False}}
                ]
            })
        else:
            query.append({"sample_sheet.group": {"$in": group}})

    if len(query) == 0:
        query = {}
    else:
        query = {"$and": query}

    if qc_list is not None and run_name is not None and len(qc_list) != 0:
        #pass
        query_result = filter_qc(db, qc_list, query)
    else:
        query_result = list(db.samples.find(query, projection)
                            .sort([(spe_field, pymongo.ASCENDING), ("name", pymongo.ASCENDING)]))
    return query_result


def get_results(sample_ids):
    connection = get_connection()
    db = connection.get_database()
    return list(db.sample_components.find({
        "sample._id": {"$in": sample_ids},
        "summary": {"$exists": True},
        "status": "Success",
        "component.name": {"$ne": "qcquickie"} #Saving transfers
    }, {"summary": 1, "sample._id": 1, "component.name" : 1, "setup_date": 1, "status": 1}).sort([("setup_date", 1)]))


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


# Run_checker.py
def get_sample_component_status(sample_ids):
    with get_connection() as connection:
        db = connection.get_database()
        sample_ids = list(map(lambda x: ObjectId(x), sample_ids))
        s_c_list = list(db.sample_components.aggregate([
            {
                "$match": {
                    "sample._id": {
                        "$in": sample_ids
                    }
                }
            },
            {
                "$sort": SON([("setup_date", 1)])
            },
            {
                "$group": {
                    "_id": {
                        "sample": "$sample._id",
                        "component": "$component.name"
                    },
                    "status": {
                        "$last": "$status"
                    }
                }
            },
            {
                "$group": {
                    "_id": "$_id.sample",
                    "s_cs": {
                        "$push": {
                            "k": "$_id.component",
                            "v": "$status"
                        }
                    }
                }
            },
            {
                "$project": {
                    "_id": 1.0,
                    "s_cs": {
                        "$arrayToObject": "$s_cs"
                    }
                }
            }
        ]))
        return s_c_list



def get_species_QC_values(ncbi_species):
    connection = get_connection()
    db = connection.get_database('bifrost_species')
    if ncbi_species != "default":
        return db.species.find_one({"ncbi_species": ncbi_species}, {"min_length": 1, "max_length": 1})
    else:
        return db.species.find_one({"organism": ncbi_species}, {"min_length": 1, "max_length": 1})


def get_sample_QC_status(last_runs):
    connection = get_connection()
    db = connection.get_database()
    samples = [sample
                for run in last_runs
                for sample in run["samples"]]
        
    samples_full = get_samples(list(map(lambda x: x["_id"], samples)))
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
                    sample_db = db.samples.find_one(
                        {"_id": run_sample["_id"]}, {"reads": 1, "stamps": 1})
                    if sample_db is not None:
                        stamps = sample_db.get("stamps", {})
                        qc_val = stamps.get(
                            "ssi_stamper", {}).get("value", "N/A")
                        if qc_val == "N/A" and (not "reads" in sample_db or not "R1" in sample_db["reads"]):
                            qc_val = "CF(LF)"
                        expert_check = False
                        if "supplying_lab_check" in stamps and "value" in stamps["supplying_lab_check"]:
                            qc_val = stamps["supplying_lab_check"]["value"]
                            expert_check = True

                        if qc_val == "fail:supplying lab":
                            qc_val = "SL"
                        elif (qc_val == "fail:core facility" or
                                qc_val == "fail:resequence"):
                            qc_val = "CF"
                        elif qc_val == "pass:OK" or qc_val == "pass:accepted":
                            qc_val = "OK"

                        if expert_check:
                            qc_val += "*"
                        sample_dict[run["name"]] = qc_val
        samples_runs_qc[name] = sample_dict
    return samples_runs_qc


def get_last_runs(run, n):
    connection = get_connection()
    db = connection.get_database()
    return list(db.runs.find({"name": {"$lte": run},
                              }, #"type": "routine"},
                             {"name": 1, "samples": 1}).sort([("_id", pymongo.DESCENDING)]).limit(n))


def get_sample(sample_id):
    connection = get_connection()
    db = connection.get_database()
    return db.samples.find_one({"_id": sample_id})


def save_sample_to_file(data_dict):
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


def get_run(run_name):
    connection = get_connection()
    db = connection.get_database()
    return db.runs.find_one({"name": run_name})


def get_sample(sample_id):
    connection = get_connection()
    db = connection.get_database()
    return db.samples.find_one({"_id": sample_id})


def get_samples(sample_ids):
    connection = get_connection()
    db = connection.get_database()
    return list(db.samples.find({"_id": {"$in": sample_ids}}))


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
