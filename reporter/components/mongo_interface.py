import pymongo
from bson.objectid import ObjectId
import keys  # .gitgnored file

def get_connection():
    "Return mongodb connection"
    return pymongo.MongoClient(keys.mongodb_url)

def get_samples_from_run_id(run_id):
    """Returns list of samples with a specific run id
    """
    with get_connection() as connection:
        db = connection.get_default_database()
        samples = db.samples
        run_samples = []
        for sample in samples.find({"run_id": run_id}):
            run_samples.append(sample)
        return run_samples


def test_get_all_samples():
    """Returns list of all samples, this function is used for testing
    """
    with get_connection() as connection:
        db = connection.get_default_database()
        samples = db.samples
        run_samples = []
        for sample in samples.find():
            run_samples.append(sample)
        return run_samples


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
        run = db.runs.find({"run.name": name}).limit(1).count(True)
    return run is not 0

def get_run_list():
    with get_connection() as connection:
        db = connection.get_default_database()
        # Fastest.
        # Change run type when this is fixed in the db
        runs = list(db.runs.find({"run.type": "project"},
                                 {"run.name": 1,
                                  "_id": 0,
                                  "samples": 1}))
    return runs

def get_group_list(run_name=None):
    with get_connection() as connection:
        db = connection.get_default_database()
        if run_name is not None:
            sample_ids = list(db.runs.find_one( #Should be taking _id
                {"run.name": run_name},
                {
                    "_id": -1,
                    "samples" : 1
                }
            )["samples"].keys())
            groups = list(db.samples.aggregate([
                {
                    "$match": {
                        "sample.name": {"$in": sample_ids},  # Should be id
                        "sample.run_folder" : {"$regex" : run_name} # Delete when we switch to id
                    }
                },
                {
                    "$group": {
                        "_id": "$sample.sample_sheet.group",
                        "count": {"$sum": 1}
                    }
                }
            ]))
        else:
            groups = list(db.samples.aggregate([
                {
                    "$group": {
                        "_id": "$sample.sample_sheet.group",
                        "count": { "$sum":1 }
                    }
                }
            ]))

    return groups


def get_species_list(run_name=None):
    with get_connection() as connection:
        db = connection.get_default_database()
        if run_name is not None:
            sample_ids = list(db.runs.find_one(  # Should be taking _id
                {"run.name": run_name},
                {
                    "_id": -1,
                    "samples": 1
                }
            )["samples"].keys())
            species = list(db.samples.aggregate([
                {
                    "$match": {
                        "sample.name": {"$in": sample_ids},  # Should be id
                        # Delete when we switch to id
                        "sample.run_folder": {"$regex": run_name}
                    }
                },
                {
                    "$group": {
                        "_id": "$qcquickie.summary.name_classified_species_1",
                        "count": {"$sum": 1}
                    }
                }
            ]))
        else:
            species = list(db.samples.aggregate([
                {
                    "$group": {
                        "_id": "$qcquickie.summary.name_classified_species_1",
                        "count": {"$sum": 1}
                    }
                }
            ]))

    return species
