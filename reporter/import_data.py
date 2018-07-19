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
