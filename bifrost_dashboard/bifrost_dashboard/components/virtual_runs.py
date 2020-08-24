import re
import datetime

import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
from flask import request
import bifrostapi



def virtual_run_options():
    runs = bifrostapi.get_run_list(run_type="virtual")
    return [{"value": x["name"], "label": x["name"]} for x in runs]

def run_full_name(name):
    name = name.upper()
    day = datetime.date.today().strftime("%y%m%d")
    return f"{day}_CUSTOM_SSIVRT_{name}"


def validate_virtual_run(samples, name, new_run=False):
    if len(name) == 0:
        return "You must enter a run name"
    if new_run:
        if re.match("^[a-zA-Z0-9_-]+$", name) is None:
            return "Invalid name. Name should only contain alphanumeric characters, - or _"
        name = run_full_name(name)
        if bifrostapi.check_run_name(name):
            return f"Run with name {name} already exists in db."
    else:
        run = bifrostapi.get_run(name)
        if run is None:
            return "Invalid run name. Run does not exist."
        elif run["type"] != "virtual":
            return "Run is not a virtual one."
    if len(samples) == 0:
        return "You must select at least one sample."
    return None

def create_run(samples, name):
    message = validate_virtual_run(samples, name, new_run=True)
    if message is not None:
        return message, False

    ip = request.remote_addr
    name = run_full_name(name)
    try:
        bifrostapi.create_virtual_run(name, ip, samples)
    except ValueError:
        return "Run could not be created", False

    return ["Run created successfully. You can access it ", html.A("here", href="/collection/" + name)]

def add_samples_to_virtual_run(samples, name):
    message = validate_virtual_run(samples, name, new_run=False)
    if message is not None:
        return message
    ip = request.remote_addr
    try:
        bifrostapi.add_samples_to_virtual_run(name, ip, samples)
    except ValueError:
        return "You can only edit virtual runs."
    return ["Samples added successfully. You can access the run ", html.A("here", href="/collection/" + name)]


def remove_samples_from_virtual_run(samples, name):
    

    ip = request.remote_addr
    try:
        error = bifrostapi.remove_samples_from_virtual_run(name, ip, samples)
    except ValueError:
        return "You can only edit virtual runs."
    if error is not None:
        return error
    return "Samples removed successfully. Refresh the page to see the changes."


def modal_add():
    return dbc.Modal(
            [
            dbc.ModalHeader("Add samples to collection"),
            dbc.ModalBody([
                dbc.Form([
                    dbc.FormGroup([
                        dbc.Label("Existing virtual run ", html_for="virtual-run-selector"),
                        dcc.Dropdown(
                                    id="virtual-run-selector",
                                    value="_new",
                                    options=[{"value": "_new", "label": "New virtual run (enter name below)"}] + virtual_run_options()
                                ),
                        dbc.Label("New virtual run ", html_for="virtual_run_name", className="mt-1"),
                        dbc.Input(type="text", id="virtual_run_name", pattern="^[a-zA-Z0-9_-]+$"),
                        dbc.FormText(["Only accepts 'A-Za-z0-9-_'"] ,color="secondary"),
                        dbc.FormText(["Full name will be: ", html.Span(id="virtual_run_full_name")] ,color="secondary")
                    ])
                ]),
                html.Div([], id="selection-modal-status")
            ]),
            dbc.ModalFooter([
                dbc.Button("Close", id="add-collection-modal-close", className="ml-auto"),
                dbc.Button("Create virtual run", id="selection-modal-create", 
                        color="primary", n_clicks=0)
            ]),
        ],
        id="add-collection-modal",
    )

def modal_remove():
    return dbc.Modal(
    [
        dbc.ModalHeader("Remove samples from collection"),
        dbc.ModalBody([
            html.P("Remove selected samples from current collection?"),
            html.Div([], id="remove-modal-status")
        ]),
        dbc.ModalFooter([
            dbc.Button("Close", id="remove-collection-modal-close", className="ml-auto"),
            dbc.Button("Remove samples", id="remove-collection-modal-remove", 
                       color="primary", n_clicks_timestamp=-1)
        ]),
    ],
    id="remove-collection-modal",
)