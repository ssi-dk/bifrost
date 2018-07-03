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
        db = connection.ngs_runs
        samples = db.samples
        run_samples = []
        for sample in samples.find({"run_id": run_id}):
            run_samples.append(sample)
        return run_samples


def test_get_all_samples():
    """Returns list of all samples, this function is used for testing
    """
    with get_connection() as connection:
        db = connection.ngs_runs
        samples = db.samples
        run_samples = []
        for sample in samples.find():
            run_samples.append(sample)
        return run_samples


def get_species_colors():
    "Get a dict with ncbi species name and color"
    with get_connection() as connection:
        db = connection.ngs_runs
        species = db.species
        colors = {}
        for species in species.find():
            colors[species["ncbi_species"]] = species["color"]
    return colors
