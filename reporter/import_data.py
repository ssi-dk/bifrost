import pandas as pd
import datetime
import components.mongo_interface as mongo_interface
from pandas.io.json import json_normalize
# Utils

def get_from_path(path_string, response):
    fields = path_string.split(".")
    for field in fields:
        response = response.get(field, None)
        if response is None:
            return response
    return response

# Main functions


def extract_data(sample_db):
    sample = {}
    sample['_id'] = str(sample_db['_id'])
    sample["input_read_status"] = sample_db["sample"]["input_read_status"]
    sample['name'] = sample_db['sample'].get('name')
    sample['user'] = sample_db['sample'].get('user')
    sample['R1_location'] = sample_db['sample'].get('R1')
    sample['R2_location'] = sample_db['sample'].get('R2')
    sample["run_name"] = sample_db["sample"]["run_folder"].split("/")[-1]
    if "setup_time" in sample_db["sample"]:
        sample["setup_time"] = datetime.datetime.strptime(
            sample_db["sample"]["setup_time"], "%Y-%m-%d %H:%M:%S.%f")
    else:
        sample["setup_time"] = None

    # Not nice, but we have many different combinations
    if "sample_sheet" in sample_db["sample"]:
        sample['supplied_name'] = sample_db['sample']["sample_sheet"]["sample_name"]
        if sample['name'] == None:
            sample['name'] = sample["supplied_name"]
        sample["provided_species"] = sample_db['sample']["sample_sheet"].get(
            "provided_species", "")
        sample['supplying_lab'] = sample_db['sample']["sample_sheet"]['group']
        sample['comments'] = sample_db["sample"]["sample_sheet"]["Comments"]
        sample["emails"] = sample_db["sample"]["sample_sheet"]["emails"]

    if "qcquickie" in sample_db:
        for key, value in sample_db['qcquickie']['summary'].items():
            sample["qcquickie_" + key] = value
        sample["qcquickie_N50"] = sample_db["qcquickie"]["quast/report_tsv"]["N50"]
        sample["qcquickie_N75"] = sample_db["qcquickie"]["quast/report_tsv"]["N75"]
        sample["qcquickie_bin_length_1x_25x_diff"] = sample["qcquickie_bin_length_at_1x"] - \
            sample["qcquickie_bin_length_at_25x"]

    if "assembly" in sample_db:
        for key, value in sample_db['assembly']['summary'].items():
            sample["assembly_" + key] = value
        if "bin_length_at_25x" not in sample_db["assembly"]["summary"]:
            print(sample["name"])
            print(sample_db['assembly']['summary'])
    
        sample["assembly_N50"] = sample_db["assembly"]["quast/report_tsv"]["N50"]
        sample["assembly_N75"] = sample_db["assembly"]["quast/report_tsv"]["N75"]
        sample["assembly_bin_length_1x_25x_diff"] = sample["assembly_bin_length_at_1x"] - \
            sample["assembly_bin_length_at_25x"]

    if not "run_name" in sample:
        sys.stderr.write("Sample {} has no run name.\n".format(sample["name"]))
        sample["run_name"] = ""

    if "qcquickie_name_classified_species_1" in sample:
        species_words = sample["qcquickie_name_classified_species_1"].split(
        )
        sample["short_class_species_1"] = '{}. {}'.format(
            species_words[0][0], ' '.join(species_words[1:]))
    else:
        sample["qcquickie_name_classified_species_1"] = "Not classified"
        sample["short_class_species_1"] = "Not classified"

    if not "supplying_lab" in sample:
        sample["supplying_lab"] = "Not specified"

    return sample

def import_data():
    return pd.DataFrame(list(map(extract_data, mongo_interface.test_get_all_samples())))


def get_species_plot_data(species_list, id_list):
    res = mongo_interface.get_species_plot_data(species_list, id_list)
    data = {}
    for doc in res:
        data[doc["_id"]] = doc["bin_coverage_at_1x"]
    return data

def check_run_name(name):
    return mongo_interface.check_run_name(name)

def get_run_list():
    return mongo_interface.get_run_list()
    
def get_group_list(run_name=None):
    return mongo_interface.get_group_list(run_name)

def get_species_list(run_name=None):
    return mongo_interface.get_species_list(run_name)


# def filter(species=None, group=None, run_name=None, aggregate=None):
#     result = mongo_interface.filter({path: 1},
#                                     run_name, species, group, aggregate)
#     return list(map(lambda x: get_from_path(path, x), result))


def filter_name(species=None, group=None, run_name=None):
    result = mongo_interface.filter({"sample.name": 1},
                                    run_name, species, group)
    return list(result)


def filter_plot(path, species=None, group=None, run_name=None, aggregate=None):

    query_result =  mongo_interface.filter(
        {
            path: 1,
            "sample.name" : 1,
            "qcquickie.summary.name_classified_species_1" : 1
        },
        run_name, species, group, aggregate)
    clean_result = []
    for item in query_result:
        clean_result.append({
            "_id": str(item["_id"]),
            'name': item["sample"]["name"],
            'value': get_from_path(path, item),
            'species': item["qcquickie"]["summary"]["name_classified_species_1"]
        })
    return pd.DataFrame(clean_result)


def filter_all(samples=[]):
    projection = {
        'sample': 1,
        'qcquickie.summary': 1,
        'assembly.summary': 1
    }
    query_result = mongo_interface.filter(projection=projection,samples=samples)
    dataframe = json_normalize(query_result)
    if "_id" in dataframe: # False for empty results.
        dataframe['_id'] = dataframe['_id'].astype(str)
    return dataframe
