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

def filter_name(species=None, group=None, run_name=None):
    result = mongo_interface.filter({"name": 1},
                                    run_name, species, group)
    return list(result)


def filter_plot(species=None, group=None, run_name=None, func=None):

    query_result =  mongo_interface.filter(
        {
            "name" : 1,
            "properties.species" : 1
        },
        run_name, species, group)
    clean_result = {}
    sample_ids = []
    for item in query_result:
        sample_ids.append(item['_id'])
        clean_result[str(item["_id"])] = {
            "_id": str(item["_id"]),
            'name': item["name"],
            'species': item["properties"]["species"]
        }
    component_result = mongo_interface.get_results(sample_ids, 'qcquickie')
    
    for item in component_result:
        s_id = str(item["sample"]["_id"])
        if 'summary' in item:
            for summary_key, summary_value in item['summary'].items():
                clean_result[s_id]['qcquickie_' + summary_key] = summary_value
        else:
            print('Missing summary', item)
        if func is not None:
            clean_result[s_id] = func(clean_result[s_id])
    return pd.DataFrame.from_dict(clean_result, orient='index')


def filter_all(samples=[]):
    projection = {
        'sample': 1,
        'qcquickie.summary': 1,
        'assembly.summary': 1
    }
    query_result = mongo_interface.filter(projection=projection,
                                          samples=samples)
    dataframe = json_normalize(query_result)
    if "_id" in dataframe: # False for empty results.
        dataframe['_id'] = dataframe['_id'].astype(str)
    return dataframe
