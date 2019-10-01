import sys
import math
import pandas as pd
from datetime import datetime
import components.mongo_interface as mongo_interface
from pandas.io.json import json_normalize
from bson.objectid import ObjectId
import components.global_vars as global_vars
import keys
from bson.json_util import dumps, loads

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

##NOTE SPLIT/SHORTEN THIS FUNCTION
def filter_all(species=None, species_source=None, group=None,
               qc_list=None, run_names=None, sample_ids=None,
               sample_names=None,
               pagination=None,
               projection=None):
    if sample_ids is None:
        query_result = mongo_interface.filter(
            run_names=run_names, species=species,
            species_source=species_source, group=group,
            qc_list=qc_list,
            sample_names=sample_names,
            pagination=pagination,
            projection=projection)
    else:
        query_result = mongo_interface.filter(
            samples=sample_ids, pagination=pagination,
            projection=projection)
    return pd.io.json.json_normalize(query_result)


def get_assemblies_paths(samples):
    return mongo_interface.get_assemblies_paths(samples)

# For run_checker
def get_sample_component_status(samples):
    return mongo_interface.get_sample_component_status(samples)


def get_species_QC_values(ncbi_species):
    return mongo_interface.get_species_QC_values(ncbi_species)

def get_sample_QC_status(run):
    return mongo_interface.get_sample_QC_status(run)


def get_last_runs(run=None, n=12, runtype=None):
    return mongo_interface.get_last_runs(run, n, runtype)

def create_feedback_s_c(user, sample, value):
    if user == "":
        raise ValueError("Missing user value")
    component = get_component(name="user_feedback", version="1.0")
    if component is None:
        component = mongo_interface.save_component(global_vars.feedback_component)
    d = datetime.utcnow()
    now = d.replace(microsecond=math.floor(d.microsecond/1000)*1000)

    if value not in ("OK", "other", "resequence"):
        raise ValueError("Feedback value {} not valid".format(value))

    if value == "OK":
        status = "pass"
    else:
        status = "fail"

    summary = {
        "stamp": {
            "name": "user_feedback",
            "display_name": "User QC feedback",
            "value": value,
            "status": status,
            "reason": ""
        }
    }

    s_c = {
        "component": {
            "_id": component["_id"],
            "name": component["name"]
        },
        "sample": {
            "_id": sample["_id"],
            "name": sample["name"]
        },
        "metadata": {
            "created_at": now,
            "updated_at": now,
        },
        "path": None,
        "status": "Success",
        "setup_date": now,
        "results": summary,
        "properties": {
            "summary": summary,
            "component": {
                "_id": component["_id"],
                "date": now,
            }
        }
        
    }

    return s_c

def add_user_feedback_to_properties(sample, s_c):
    """
    Adds/updates user feedback to sample properties
    """
    stamper = sample.get("properties", {}).get("stamper", {})
    stamper["component"] = s_c["properties"]["component"]
    summary = stamper.get("summary", {})

    old_value = summary.get("stamp").get("value")

    summary["stamp"] = s_c["properties"]["summary"]["stamp"]
    sample["properties"] = sample.get("properties", {})
    sample["properties"]["stamper"] = stamper
    return sample, old_value


def add_user_feedback(user, sample_id, value):
    """
    Submits user feedback to a sample
    """
    sample = mongo_interface.get_samples([ObjectId(sample_id)])[0]
    s_c = create_feedback_s_c(user, sample, value)
    sample, old_value = add_user_feedback_to_properties(sample, s_c)
    mongo_interface.save_sample_component(s_c)
    mongo_interface.save_sample(sample)
    return (sample["name"], old_value)

def add_batch_user_feedback_and_mail(feedback_pairs, user):
    """
    Main function called to send sample QC feedback 
    """
    email_pairs = []
    for sample_id, value in feedback_pairs:
        name, old_value = add_user_feedback(user, sample_id, value)
        runs = mongo_interface.get_sample_runs([ObjectId(sample_id)])
        run_names = [run["name"] for run in runs]
        if "fail:core facility" in (old_value, value):
            email_pairs.append(name, old_value, value, run_names)
    send_mail(email_pairs, user)


def get_component(name, version):
    return mongo_interface.get_component(name, version)


def get_run(run_name):
    return mongo_interface.get_run(run_name)


def send_mail(sample_info, user):
    """
    Sends email about sample feedback updates.
    """
    import smtplib
    from email.mime.multipart import MIMEMultipart
    from email.mime.text import MIMEText

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
