import pymongo
import keys  # .gitgnored file
from bson.objectid import ObjectId


def get_connection():
    "Return mongodb connection"
    return pymongo.MongoClient(keys.mongodb_url)

def check_run_name(name):
    with get_connection() as connection:
        db = connection.get_database()
        # Fastest.
        run = db.runs.find({"name": name}).limit(1).count(True)
    return run is not 0

def get_run_list():
    with get_connection() as connection:
        db = connection.get_database()
        # Fastest.
        runs = list(db.runs.find({"type": "routine"}, #Leave in routine
                                 {"name": 1,
                                  "_id": 0,
                                  "samples": 1}).sort([['name', pymongo.DESCENDING]]))
    return runs

def get_group_list(run_name=None):
    """
    Get group list but most importantly count of samples per group for a run.
    """
    with get_connection() as connection:
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
    with get_connection() as connection:
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

    result = db.samples.aggregate([
        {
            "$match": {"$and" :query},
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
                                    { "$eq" : ["$component.name", "ssi_stamper"]},
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
            "$match": {"$or": qc_query}
            
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
    with get_connection() as connection:
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

        if qc_list is not None and run_name is not None and len(qc_list) != 0:
            #pass
            query_result = filter_qc(db, qc_list, query)
        else:
            query_result = list(db.samples.find({"$and": query}, projection)
                                .sort([(spe_field, pymongo.ASCENDING), ("name", pymongo.ASCENDING)]))
        return query_result



def get_results(sample_ids):
    with get_connection() as connection:
        db = connection.get_database()
        return list(db.sample_components.find({
            "sample._id": {"$in": sample_ids},
            "summary": {"$exists": True},
            "status": "Success",
            "component.name": {"$ne": "qcquickie"} #Saving transfers
        }, {"summary": 1, "sample._id": 1, "component.name" : 1, "setup_date": 1, "status": 1}).sort([("setup_date", 1)]))

def get_sample_runs(sample_ids):
    with get_connection() as connection:
        db = connection.get_database()
        return list(db.runs.find({"samples": {"$elemMatch": {"_id": {"$in": sample_ids}}}}))


def get_assemblies_paths(sample_ids):
    with get_connection() as connection:
        db = connection.get_database()
        return list(db.sample_components.find({
            "sample._id": {"$in": list(map(lambda x: ObjectId(x), sample_ids))},
            "component.name": "assemblatron"
        }, {"path": 1, "sample._id": 1}))

# Run_checker.py
def get_sample_component_status(run_name):
    with get_connection() as connection:
        db = connection.get_database()
        run = db.runs.find_one({"name": run_name})
        samples_ids = list(map(lambda x:x["_id"], run["samples"]))
        components_names = list(map(lambda x: x["name"], run["components"]))
        s_c_list = db.sample_components.find({
            "sample._id": {"$in": samples_ids},
            "component.name": {"$in": components_names}
        }, {"sample": 1, "status": 1, "component.name": 1})
        output = {}
        for s_c in s_c_list:
            sample = output.get(s_c["sample"]["name"], {
                "sample._id": s_c["sample"]["_id"]
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
            output[s_c["sample"]["name"]] = sample
        return output


def get_species_QC_values(ncbi_species):
    with get_connection() as connection:
        db = connection.get_database()
        if ncbi_species != "default":
            return db.species.find_one({"ncbi_species": ncbi_species}, {"min_length": 1, "max_length": 1})
        else:
            return db.species.find_one({"organism": ncbi_species}, {"min_length": 1, "max_length": 1})


def get_sample_QC_status(last_runs):
    with get_connection() as connection:
        db = connection.get_database()
        samples = [sample
                   for run in last_runs
                for sample in run["samples"]]

        samples_runs_qc = {}
        for sample in samples:
            sample_dict = {}

            name = sample["name"]
            for run in last_runs:
                for run_sample in run["samples"]:
                    if name == run_sample["name"]:
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
                            elif qc_val == "fail:core facility":
                                qc_val = "CF"
                            elif qc_val == "pass:OK":
                                qc_val = "OK"

                            if expert_check:
                                qc_val += "*"
                            sample_dict[run["name"]] = qc_val
            samples_runs_qc[sample["name"]] = sample_dict
        return samples_runs_qc

def get_last_runs(run, n): #merge with  get_run_list
    with get_connection() as connection:
        db = connection.get_database()
        return list(db.runs.find({"name": {"$lte": run}, "type": "routine"}, {"name": 1, "samples": 1}).sort([("name", pymongo.DESCENDING)]).limit(n))


def get_samples(sample_ids):
    sample_ids = list(map(lambda x:ObjectId(x), sample_ids))
    with get_connection() as connection:
        db = connection.get_database()
        return list(db.samples.find({"_id": {"$in": sample_ids}}))
        

def save_sample(data_dict):
    """COPIED FROM BIFROSTLIB Insert sample dict into mongodb.
    Return the dict with an _id element"""
    with get_connection() as connection:
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

def get_run(run_name): #to be merged with get_run_list
    with get_connection() as connection:
        db = connection.get_database()
        return db.runs.find_one({"name": run_name})
