import os
import pymongo
import keys  # .gitgnored file
from bson.objectid import ObjectId
import atexit


PAGESIZE = 25

CONNECTION = None
SPECIES_CONNECTION = None

TABLE_QUERY = [
    {
        '$lookup': {
            'from': 'sample_components',
            'let': {
                'sample_id': '$_id'
            },
            'pipeline': [
                {
                    '$match': {
                        'status': 'Success',
                        'name': {
                            '$nin': [
                                'qcquickie', 'analyzer'
                            ]
                        },
                        '$expr': {
                            '$eq': [
                                '$sample._id', '$$sample_id'
                            ]
                        }
                    }
                }, {
                    '$project': {
                        'component.name': 1,
                        'setup_date': 1,
                        'summary': 1,
                        'status': 1,
                        '_id': 0
                    }
                }, {
                    '$sort': {
                        'component.name': 1,
                        'setup_date': -1
                    }
                }, {
                    '$group': {
                        '_id': '$component.name',
                        's_c_g': {
                            '$push': '$$ROOT'
                        }
                    }
                }, {
                    '$project': {
                        '_id': 0,
                        'v': {
                            '$arrayElemAt': [
                                '$s_c_g', 0
                            ]
                        },
                        'k': '$_id'
                    }
                }
            ],
            'as': 'sample_components'
        }
    }, {
        '$project': {
            'name': 1,
            '_id': 1,
            'component.name': 1,
            'sample_sheet': 1,
            'stamps': 1,
            'properties': 1,
            'sample_components': {
                '$arrayToObject': '$sample_components'
            }
        }
    }
]


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
    runs = list(db.runs.find({"type": "routine"}, #Leave in routine
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


def get_species_list(species_source):
    connection = get_connection()
    db = connection.get_database()
    if species_source == "provided":
        spe_field = "properties.provided_species"
    else:
        spe_field = "properties.detected_species"
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
        runs = list(db.runs.find({"type": "routine"}, {"samples": 1}))
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


def filter_qc(qc_list):
    if qc_list is None or len(qc_list) == 0:
        return None
    qc_query = []
    for elem in qc_list:
        if elem == "Not checked":
            qc_query.append({"$and": [
                {"reads.R1": {"$exists": True}},
                {"stamps.ssi_stamper": {"$exists": False}}
            ]})
        elif elem == "fail:core facility":
            qc_query.append({"$or": [
                        {"reads.R1": {"$exists": False}},
                        {"stamps.ssi_stamper.value": "fail:core facility"}
                    ]
                })
        else:
            qc_query.append({"stamps.ssi_stamper.value": elem})
    return {"$match": {"$and": qc_query}}


def filter(projection=None, run_names=None,
           species=None, species_source="species", group=None,
           qc_list=None, samples=None, pagination=None,
           sample_names=None):
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
    if sample_names is not None and len(sample_names) != 0:
        query.append({"name": {"$in": sample_names}})
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
                    {"sample_sheet.group": None},
                    {"sample_sheet.group": {"$in": group}},
                    {"sample_sheet.group": {"$exists": False}}
                ]
            })
        else:
            query.append({"sample_sheet.group": {"$in": group}})

    sort_step = {"$sort": {"name": 1}}

    if pagination is not None:
        p_limit = pagination['page_size']
        p_skip = pagination['page_size'] * pagination['current_page']
    else:
        p_limit = 1000
        p_skip = 0

    skip_limit_steps = [{"$skip": p_skip}, {"$limit": p_limit}]

    qc_query = filter_qc(qc_list)

    if len(query) == 0:
        if qc_query is None:
            final_query = [sort_step] + TABLE_QUERY + skip_limit_steps
        else:
            final_query = TABLE_QUERY + [qc_query, sort_step] + skip_limit_steps
    else:
        if qc_query is None:
            final_query = [{"$match": {"$and": query}}] + \
                TABLE_QUERY + [sort_step] + skip_limit_steps
        else:
            final_query = [{"$match": {"$and": query}}] + \
                TABLE_QUERY + [qc_query, sort_step] + skip_limit_steps



    query_result = db.samples.aggregate(final_query)
    #query_count = query_result.count()  # Count ignores limit
    query_count = 100
    query_result = list(query_result)
    return (query_result, query_count)


def get_results(sample_ids):
    connection = get_connection()
    db = connection.get_database()
    return list(db.sample_components.find({
        "sample._id": {"$in": sample_ids},
        "summary": {"$exists": True},
        "status": "Success",
        "component.name": {"$nin": ["qcquickie", "testomatic"]} #Saving transfers
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
    }, {"path": 1, "sample.name": 1}))


# Run_checker.py
def get_sample_component_status(sample_ids):
    with get_connection() as connection:
        db = connection.get_database()
        sample_ids = list(map(lambda x: ObjectId(x), sample_ids))
        s_c_list = db.sample_components.find({
            "sample._id": {"$in": sample_ids},
        }, {"sample._id": 1, "status": 1, "component.name": 1})
        output = {}
        for s_c in s_c_list:
            sample = output.get(str(s_c["sample"]["_id"]), {
                "sample._id": str(s_c["sample"]["_id"])
            })
            status = s_c["status"]
            if status == "Success":
                status = "OK"
                status_code = 2
            elif status == "Running":
                status_code = 1
            elif status == "initialized":
                status = "init."
                status_code = 0
            elif status == "Failure":
                status = "Fail"
                status_code = -1
            elif status == 'Requirements not met':
                status = "Req."
                status_code = -2
            elif status == 'queued to run':
                status = "queue"
                status_code = 0
            else:
                status_code = float('nan')
            sample[s_c["component"]["name"]] = (status_code, status)
            output[str(s_c["sample"]["_id"])] = sample
        return output


def get_species_QC_values(ncbi_species):
    connection = get_species_connection()
    db = connection.get_database()
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
        
    samples_full = get_samples(list(map(lambda x: str(x["_id"]), samples)))
    samples_by_ids = {str(s["_id"]): s for s in samples_full}

    samples_runs_qc = {}
    for sample in samples:
        sample_dict = {}

        name = samples_by_ids[str(sample["_id"])]["name"]
        for run in last_runs:
            for run_sample in run["samples"]:
                if name == samples_by_ids[str(run_sample["_id"])]["name"]:
                    stamps = db.samples.find_one(
                        {"_id": run_sample["_id"]}, {"stamps": 1})
                    if stamps is not None:
                        stamps = stamps.get("stamps", {})
                        qc_val = stamps.get(
                            "ssi_stamper", {}).get("value", "N/A")
                        expert_check = False
                        if "ssi_expert_check" in stamps and "value" in stamps["ssi_expert_check"]:
                            qc_val = stamps["ssi_expert_check"]["value"]
                            expert_check = True

                    if qc_val == "fail:supplying lab":
                        qc_val = "SL"
                    elif (qc_val == "fail:core facility" or
                            qc_val == "fail:resequence"):
                        qc_val = "CF"
                    elif qc_val == "pass:OK":
                        qc_val = "OK"

                        if expert_check:
                            qc_val += "*"
                        sample_dict[run["name"]] = qc_val
        samples_runs_qc[name] = sample_dict
    return samples_runs_qc


def get_last_runs(run, n):
    connection = get_connection()
    db = connection.get_database()
    return list(db.runs.find({"name": {"$lte": run}, "type": "routine"}, {"name": 1, "samples": 1}).sort([("name", pymongo.DESCENDING)]).limit(n))


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
