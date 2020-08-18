import re
import datetime

import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_table
from flask import request
import bifrostapi

modal = dbc.Modal(
    [
        dbc.ModalHeader("Selected samples"),
        dbc.ModalBody([
            dash_table.DataTable(id="selected-samples-table",
                                data=[],
                                columns=[{"id": "id", "name": "id"},{"id": "name", "name": "name"}],
                                row_deletable=True),
            html.Hr(),
            dbc.Form([
                dbc.FormGroup([
                    dbc.Label("Virtual run name ", html_for="virtual_run_name"),
                    dbc.Input(type="text", id="virtual_run_name", pattern="^[a-zA-Z0-9_-]+$"),
                    dbc.FormText(["Only accepts 'A-Za-z0-9-_'"] ,color="secondary"),
                    dbc.FormText(["Full name will be: ", html.Span(id="virtual_run_full_name")] ,color="secondary")
                ])
            ]),
            html.Div([], id="selection-modal-status")
        ]),
        dbc.ModalFooter([
            dbc.Button("Close", id="selection-modal-close", className="ml-auto"),
            dbc.Button("Create virtual run", id="selection-modal-create", 
                       color="primary", n_clicks_timestamp=-1, disabled=True)
        ]),
    ],
    id="selection-modal",
)

def create_run(samples, name):
    if len(name) == 0:
        return "You must enter a run name", False
    name = name.upper()
    if re.match("^[a-zA-Z0-9_-]+$", name) is None:
        return "Invalid name. Name should only contain alphanumeric characters, - or _", False
    if len(samples) == 0:
        return "Run must contain at least one sample", False

    day = datetime.date.today().strftime("%y%m%d")
    name = f"{day}_CUSTOM_SSIVRT_{name}"

    ip = request.remote_addr
    try:
        bifrostapi.create_virtual_run(name, ip, samples)
    except ValueError:
        return "Run could not be created, check that the name is not taken.", False

    return [["Run created successfully. You can access it ", html.A("here", href="/collection/" + name)], True]