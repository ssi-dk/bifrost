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
    print('group: ', group)
    result = mongo_interface.filter({"name": 1},
                                    run_name, species, group)
    return list(result)

##NOTE SPLIT/SHORTEN THIS FUNCTION
def filter_all(species=None, group=None, run_name=None, func=None, sample_ids=None):
    if sample_ids is None:
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
    else:
        query_result = mongo_interface.filter(
            {
                "name": 1,
                "properties.species": 1,
                "sample_sheet" : 1
            },
            samples=sample_ids)
        clean_result = {}
        sample_ids = []
        for item in query_result:
            item_id = str(item["_id"])
            sample_ids.append(item['_id'])
            clean_result[item_id] = {
                "_id": item_id,
                'name': item["name"],
                'species': item["properties"]["species"]
            }
            if 'sample_sheet' in item:
                for sheet_key, sheet_value in item['sample_sheet'].items():
                    clean_result[item_id]['sample_sheet.' +
                                       sheet_key] = sheet_value

    component_result = mongo_interface.get_results(sample_ids)
    for item in component_result:
        item_id = str(item["sample"]["_id"])
        component = item['component']['name']
        if 'summary' in item:
            for summary_key, summary_value in item['summary'].items():
                clean_result[item_id][component + '.' +
                                      summary_key] = summary_value
        else:
            print('Missing summary', item)
        if func is not None:
            clean_result[item_id] = func(clean_result[item_id])
    return pd.DataFrame.from_dict(clean_result, orient='index')