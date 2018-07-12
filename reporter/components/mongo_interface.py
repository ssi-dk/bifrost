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
        samples = db.samples_test
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
        samples = db.samples_test
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
        # for data in res:
        #     plot_data[str(data["_id"])] = {
        #         "assembly": {
        #             "contig_coverage": data["assembly"]["contigs_sum_cov"]["contig_depth"]
        #         },
        #         "qcquickie": {
        #             "contig_coverage": data["qcquickie"]["contigs_sum_cov"]["contig_depth"]
        #         }
        #     }
    return res

