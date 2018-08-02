import pymongo
import keys  # .gitgnored file
from bson.objectid import ObjectId
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
        db = connection.get_default_database()
        species_col = db.species
        colors = {}
        for species in species_col.find():
            colors[species["organism"]] = species["color"]
    return colors

def get_species_plot_data(species_list, id_list):
    """Get plot data for many samples using a list of Ids"""
    with get_connection() as connection:
        db = connection.get_default_database()
        samples = db.samples
        plot_data = {}
        id_list = list(map(lambda x: ObjectId(x), id_list))
        res = db.samples.aggregate([
            {"$match": {
                "qcquickie.summary.name_classified_species_1": {"$in": species_list},
                "_id": {"$not": {"$in": id_list}}
                }
            },
            {"$group": {
                "_id": "$qcquickie.summary.name_classified_species_1",
                "bin_coverage_at_1x": {"$push": "$qcquickie.summary.bin_length_at_1x"}
                }
            }
        ])
    return res

def check_run_name(name):
    with get_connection() as connection:
        db = connection.get_default_database()
        # Fastest.
        run = db.runs.find({"name": name}).limit(1).count(True)
    return run is not 0

def get_run_list():
    with get_connection() as connection:
        db = connection.get_default_database()
        # Fastest.
        runs = list(db.runs.find({"type": "testing"}, #Leave in routine
                                 {"name": 1,
                                  "_id": 0,
                                  "samples": 1}))
    return runs

def get_group_list(run_name=None):
    with get_connection() as connection:
        db = connection.get_default_database()
        if run_name is not None:
            run_samples = db.runs.find_one(
                {"name": run_name},
                {
                    "_id": 0,
                    "samples._id": 1
                }
            )['samples']
            sample_ids = [s['_id'] for s in run_samples]
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
        db = connection.get_default_database()
        if run_name is not None:
            run_samples = db.runs.find_one(
                {"name": run_name},
                {
                    "_id": 0,
                    "samples._id": 1
                }
            )['samples']
            sample_ids = [s['_id'] for s in run_samples]
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
                }
            ]))
        else:
            species = list(db.samples.aggregate([
                {
                    "$group": {
                        "_id": "$properties.species",
                        "count": {"$sum": 1}
                    }
                }
            ]))

    return species


def filter(projection=None, run_name=None,
           species=None, group=None, samples=None):
    with get_connection() as connection:
        db = connection.get_default_database()
        query = []
        sample_set = set()
        if samples is not None:
            sample_set = {ObjectId(id) for id in samples}
            query.append({"_id": {"$in": list(sample_set)}})
        if run_name is not None and run_name != "":
            run_samples = db.runs.find_one(
                {"name": run_name},
                {
                    "_id": 0,
                    "samples._id": 1
                }
            )['samples']
            run_sample_set = {s['_id'] for s in run_samples}
        
            if len(sample_set):
                inter = run_sample_set.intersect(sample_set)
                query.append({"_id": {"$in": list(inter)}})
            else:
                query.append({"_id": {"$in": list(run_sample_set)}})
        if species is not None and len(species) != 0:
            query.append({"$or":
                [
                    {"properties.species": None},
                    {"properties.species": {"$in": species}},
                    {'properties.species': {'$exists': False}}
                ]
            })
        if group is not None and len(group) != 0:
            if "Not defined" in group:
                query.append({"$or":
                    [
                        {"sample_sheet.group": None},
                        {"sample_sheet.group": {"$in": group}},
                        {'sample_sheet.group': {'$exists': False}}
                    ]
                })
            else:
                query.append({"sample_sheet.group": {"$in": group}})
        return list(db.samples.find({'$and': query}, projection))


def get_results(sample_ids):
    with get_connection() as connection:
        db = connection.get_default_database()
        return list(db.sample_components.find({
            'sample._id': {'$in': sample_ids}
        }, {'summary': 1, 'sample._id': 1, 'component.name' : 1}))

