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

def get_plot_data(id_list):
    """Get plot data for many samples using a list of Ids"""
    with get_connection() as connection:
        db = connection.get_default_database()
        samples = db.samples_test
        plot_data = {}
        res = samples.find(
            {
                "_id": {
                    "$in": [ObjectId(id) for id in id_list]
                }
            }
        )
        for data in res:
            plot_data[str(data["_id"])] = {
                "assembly": {
                    "contig_coverage": data["assembly"]["contigs_sum_cov"]["contig_depth"]
                },
                "qcquickie": {
                    "contig_coverage": data["qcquickie"]["contigs_sum_cov"]["contig_depth"]
                }
            }
    return plot_data

def _format_contigs(contig_dict):
    contig_list = []
    for key, value in contig_dict.items():
        contig_list.append({
            "name": key,
            "coverage": value["coverage"],
            "total_depth": value["total_depth"],
            "total_length": value["total_length"]
        })
    return contig_list
