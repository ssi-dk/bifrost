import sys
import pandas as pd
from datetime import datetime
import components.mongo_interface as mongo_interface
from pandas.io.json import json_normalize
from bson.objectid import ObjectId
import components.global_vars as global_vars
import keys

pd.options.mode.chained_assignment = None


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

def get_species_list(species_source, run_name=None):
    return mongo_interface.get_species_list(species_source, run_name)

def filter_name(species=None, group=None, qc_list=None, run_name=None):
    result = mongo_interface.filter({"name": 1, "sample_sheet.sample_name": 1},
                                    run_name, species, group, qc_list)
    return list(result)

def filter_all(species=None, species_source=None, group=None, qc_list=None, run_name=None, func=None, sample_ids=None):

    if sample_ids is None:
        query_result =  mongo_interface.filter(
            {
                "name" : 1,
                "properties": 1,
                "sample_sheet": 1,
                "report": 1,
                "reads": 1,
                "stamps": 1
            },
            run_name, species, species_source, group, qc_list=qc_list)
    else:
        query_result = mongo_interface.filter(
            {
                "name": 1,
                "properties": 1,
                "sample_sheet": 1,
                "report": 1,
                "reads": 1,
                "stamps": 1
            },
            samples=sample_ids)

    # Check if any sample belongs to the old schema.
    old_schema = []
    new_schema = []
    for result in query_result:
        if not result.get("report"):
            old_schema.append(result)
        else:
            new_schema.append(result)
    if not old_schema:
        return query_result
    # clean_result = {}
    sample_ids = [x["_id"] for x in old_schema]
    samples_by_ids = {str(x["_id"]):x for x in old_schema}
    # unnamed_count = 0

    category_dict = {
        "assemblatron": "denovo_assembly",
        "ariba_mlst": "mlst",
        "ariba_resfinder": "resistance",
        "cge_mlst": "mlst",
        "cge_resfinder": "resistance",
        "min_read_check": "size_check",
        "ssi_stamper": "stamper",
        "whats_my_species": "species_detection"
    }

    component_result = mongo_interface.get_results(sample_ids)
    for item in component_result:
        item_id = str(item["sample"]["_id"])
        component = item["component"]["name"]
        if component in category_dict:
            category = category_dict[component]
            if "summary" in item:
                samples_by_ids[item_id]["properties"][category] = {
                    "summary": item["summary"]
                }
    for sample in samples_by_ids:
        new_schema.append(samples_by_ids[sample])

    return new_schema



def add_sample_runs(sample_df):
    """Returns the runs each sample belongs to"""
    sample_ids = sample_df["_id"].tolist()
    sample_ids = list(map(lambda x: ObjectId(x), sample_ids))
    runs = mongo_interface.get_sample_runs(sample_ids)
    sample_runs = {}
    # Weekend challenge: turn this into a double nested dictionary comprehension
    for run in runs:
        for sample in run["samples"]:
            if sample["_id"] in sample_ids:
                s = sample_runs.get(str(sample["_id"]), [])
                s.append(run["name"])
                sample_runs[str(sample["_id"])] = s
    sample_df.loc[:, 'runs'] = sample_df["_id"].map(sample_runs)
    return sample_df

def get_assemblies_paths(samples):
    return mongo_interface.get_assemblies_paths(samples)

# For run_checker
def get_sample_component_status(sample_ids):
    return mongo_interface.get_sample_component_status(sample_ids)


def get_species_QC_values(ncbi_species):
    return mongo_interface.get_species_QC_values(ncbi_species)

def get_sample_QC_status(run):
    return mongo_interface.get_sample_QC_status(run)

def get_last_runs(run, n):
    return mongo_interface.get_last_runs(run, n)


def post_stamps(stamplist):
    for pair in stamplist:
        sample_id, stamp = pair
        sample_db = mongo_interface.get_samples([ObjectId(sample_id)])[0]
        stamps = sample_db.get("stamps", {})
        stamp_list = stamps.get("stamp_list", [])
        stamp_list.append(stamp)
        stamps["stamp_list"] = stamp_list
        stamps[stamp["name"]] = stamp
        sample_db["stamps"] = stamps
        mongo_interface.save_sample_to_file(sample_db)
    return "Feedback saved"

def get_run(run_name):
    return mongo_interface.get_run(run_name)

def get_samples(sample_ids):
    return mongo_interface.get_samples(sample_ids)


def email_stamps(stamplist):
    import smtplib
    from email.mime.multipart import MIMEMultipart
    from email.mime.text import MIMEText

    user = ""
    sample_info = []
    for stamp_pair in stamplist:
        sample_id, stamp = stamp_pair
        user = stamp["user"]
        new_status = stamp["value"]
        sample = mongo_interface.get_samples([ObjectId(sample_id)])[0]
        if sample is not None:
            sample_name = sample["name"]
        else:
            sample_name = "DBERROR, ID: {}".format(sample_id)
        runs = mongo_interface.get_sample_runs([ObjectId(sample_id)])
        run_names = [run["name"] for run in runs]
        old_status = "none"
        if "stamps" in sample:
            if "supplying_lab_check" in sample["stamps"]:
                old_status = sample["stamps"]["supplying_lab_check"]["value"]
            elif "ssi_stamper" in sample["stamps"]:
                old_status = sample["stamps"]["ssi_stamper"]["value"]
        if "fail:resequence" in (old_status, new_status):
            sample_info.append((sample_name, old_status, new_status, run_names))

    if len(sample_info) == 0:
        return

    short_samples = ",".join([pair[0] for pair in sample_info])[
        :60]  # Trimmed to 60 chars
    msg = MIMEMultipart("alternative")
    msg["From"] = keys.email_from
    msg['Subject'] = 'Sample status change: "{}"'.format(short_samples)
    msg['To'] = keys.email_to

    email_text = 'Automatic message:\nUser "{}" has changed the status of the following samples:\n\nSample name                Old status            New status            Run name\n'.format(
        user)
    email_html = '<html><body>Automatic message:\nUser "{}" has changed the status of the following samples:\n\n<pre>Sample name                Old status            New status            Run name\n'.format(
        user)
    table = ""
    for pair in sample_info:
        table += "{:27s}{:22s}{:22s}{}\n".format(
            pair[0], pair[1], pair[2], ",".join(pair[3]))
    email_text += table
    email_html += table + "</pre></body></html>"
    print(email_text)
    msg.attach(MIMEText(email_text, 'plain'))
    msg.attach(MIMEText(email_html, 'html'))
    s = smtplib.SMTP('localhost')
    s.sendmail(msg["From"], msg["To"], msg.as_string())


def get_samples(sample_ids):
    sample_ids = [ObjectId(id) for id in sample_ids]
    return mongo_interface.get_samples(sample_ids)


def get_comment(run_id):
    returned = mongo_interface.get_comment(run_id)
    if returned:
        return returned.get("Comments", None)
    else:
        return returned

def set_comment(run_id, comment):
    return mongo_interface.set_comment(run_id, comment)


def load_file_from_db(file_id=None):
    return mongo_interface.load_file_from_db(file_id)
