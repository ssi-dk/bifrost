import pymongo
import keys  # .gitgnored file
from bson.objectid import ObjectId


PAGESIZE = 25

def get_from_path(path_string, response):
    fields = path_string.split(".")

    for field in fields:
        response = response.get(field, None)
        if response is None:
            return response
    return response

def flatten(l, path):
    """Flattens a list of dicts into a list"""
    return [get_from_path(path, el) for el in l]


def get_connection():
    "Return mongodb connection"
    return pymongo.MongoClient(keys.mongodb_url)

def get_species_colors(): 
    """Get a dict with ncbi species name and color"""
    with get_connection() as connection:
        db = connection.get_database()
        species_col = db.species
        colors = {}
        for species in species_col.find():
            colors[species["organism"]] = species["color"]
    return colors

def get_species_plot_data(species_list, id_list):
    """Get plot data for many samples using a list of Ids"""
    with get_connection() as connection:
        db = connection.get_database()
        samples = db.samples
        plot_data = {}
        id_list = list(map(lambda x: ObjectId(x), id_list))
        res = db.sample_components.aggregate([
            {"$match": {
                "summary.name_classified_species_1": {"$in": species_list},
                "sample._id": {"$not": {"$in": id_list}}
                }
            },
            {"$group": {
                "_id": "$summary.name_classified_species_1",
                "bin_coverage_at_1x": {"$push": "$summary.bin_length_at_1x"}
                }
            }
        ])
    return res

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


def get_species_list(run_name=None):
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
            species = list(db.samples.aggregate([
                {
                    "$match": {
                        "_id": {"$in": sample_ids}
                    }
                },
                {
                    "$group": {
                        "_id": "$properties.species",
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
                        "_id": "$properties.species",
                        "count": {"$sum": 1}
                    }
                },
                {
                    "$sort": {"_id": 1}
                }
            ]))

    return species


def get_qc_list(run_name=None):
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
            qcs = list(db.sample_components.aggregate([
                {
                    "$match": {
                        "sample._id": {"$in": sample_ids}
                    }
                },
                {
                    "$group": {
                        "_id": "$summary.assemblatron:action",
                        "count": {"$sum": 1}
                    }
                },
                {
                    "$sort": {"_id": 1}
                }
            ]))
        else:
            qcs = list(db.samples.aggregate([
                {
                    "$group": {
                        "_id": "$summary.assemblatron:action",
                        "count": {"$sum": 1}
                    }
                },
                {
                    "$sort": {"_id": 1}
                }
            ]))

    return qcs

def filter_qc(db, qc_list, results):
    sample_ids = list(map(lambda x: x["_id"], results))
    if "Not determined" in qc_list:
        qc_query = {
            "sample._id": {"$in": sample_ids},
            "$or":
                [
                    {"summary.assemblatron:action": None},
                    {"summary.assemblatron:action": {"$in": qc_list}},
                    {"summary.assemblatron:action": {"$exists": False}}
                ]
            }
    else:
        qc_query = {
            "sample._id": {"$in": sample_ids},
            "summary.assemblatron:action": {"$in": qc_list}
        }
    qc_query_res = list(db.sample_components.find(qc_query))
    qc_ids = list(map(lambda x: x["sample"]["_id"], qc_query_res))
    filtered = [result for result in results if result["_id"] in qc_ids]
    if "Not tested" in qc_list:
        not_tested_rev = db.sample_components.find({
            "sample._id": {"$in": sample_ids},
            "component.name": "testomatic"
        })
        not_tested_rev_ids = list(map(lambda x:x["sample"]["_id"], not_tested_rev))
        for sample in results:
            if sample["_id"] not in not_tested_rev_ids:
                filtered.append(sample)
    return filtered


def filter(projection=None, run_name=None,
           species=None, group=None, qc_list=None, samples=None, page=None):
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
                        {"properties.species": None},
                        {"properties.species": {"$in": species}},
                        {"properties.species": {"$exists": False}}
                    ]
                })
            else:
                query.append({"properties.species": {"$in": species}})
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

        query_result = list(db.samples.find({"$and": query}, projection)\
            .sort([("properties.species", pymongo.ASCENDING), ("name", pymongo.ASCENDING)]))

        if qc_list is not None and len(qc_list) != 0:
            query_result = filter_qc(db, qc_list, query_result)

        if page is None:
            return query_result
        else:
            skips = PAGESIZE * (page)
            return query_result[skips:skips+PAGESIZE]



def get_results(sample_ids):
    with get_connection() as connection:
        db = connection.get_database()
        return list(db.sample_components.find({
            "sample._id": {"$in": sample_ids}
        }, {"summary": 1, "sample._id": 1, "component.name" : 1}))

def get_sample_runs(sample_ids):
    with get_connection() as connection:
        db = connection.get_database()
        return list(db.runs.find({"samples": {"$elemMatch": {"_id": {"$in": sample_ids}}}}))


def get_read_paths(sample_ids):
    with get_connection() as connection:
        db = connection.get_database()
        return list(db.samples.find({"_id": {"$in": list(map(lambda x:ObjectId(x), sample_ids))}}, {"reads": 1, "name": 1}))

# Run_checker.py
def get_sample_component_status(run_name):
    with get_connection() as connection:
        db = connection.get_database()
        run = db.runs.find_one({"name": run_name})
        samples_ids = list(map(lambda x:x["_id"], run["samples"]))
        components_ids = list(map(lambda x: x["_id"], run["components"]))
        s_c_list = db.sample_components.find({
            "sample._id": {"$in": samples_ids},
            "component._id": {"$in": components_ids}
        })
        output = {}
        for s_c in s_c_list:
            sample = output.get(s_c["sample"]["name"], {})
            status = s_c["status"]
            if status == "Success":
                status_code = 2
            elif status == "initialized":
                status_code = 1
            elif status == "Failure":
                status_code = -1
            elif status == 'queued to run':
                status_code = 0
            else:
                status_code = float('nan')
            sample[s_c["component"]["name"]] = (status_code, status)
            output[s_c["sample"]["name"]] = sample
        return output

