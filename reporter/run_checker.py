# -*- coding: utf-8 -*-
import os
import sys

import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State

import dash_scroll_up

import components.mongo_interface
import import_data


#Globals
PAGESIZE = 100

app = dash.Dash()
app.config["suppress_callback_exceptions"] = True
app.title = "Serum QC Run Checker"

# Temp css to make it look nice
# Dash CSS
app.css.append_css(
    {"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})
# Lato font
app.css.append_css(
    {"external_url": "https://fonts.googleapis.com/css?family=Lato"})

app.layout = html.Div([
    html.Div(className="container", children=[
        dash_scroll_up.DashScrollUp(
            id="input",
            label="UP",
            className="button button-primary no-print"
        ),
        html.H1("SerumQC Run Checker"),
        html.H2("", id="run-name"),
        dcc.Location(id="url", refresh=False),
        html.Div(id="run-table"),
        html.Div(id="run-selector")

    ]),
    html.Footer([
        "Created with 🔬 at SSI. Bacteria icons from ",
        html.A("Flaticon", href="https://www.flaticon.com/"),
        "."], className="footer container")
], className="appcontainer")


# Callbacks

# We could make this one much faster by hiding the unused species with CSS
# by adding a new hidden class.


@app.callback(
    Output("run-name", "children"),
    [Input("url", "pathname")]
)
def update_run_name(pathname):
    if pathname is None:
        pathname = "/"
    path = pathname.split("/")
    if import_data.check_run_name(path[1]):
        return path[1]
    elif path[1] == "":
        return ""
    else:
        return "Not found"


@app.callback(
    Output("run-selector", "children"),
    [Input("run-name", "children")]
)
def update_run_list(run_name):
    if len(run_name) == 0:
        run_list = import_data.get_run_list()
        run_list_options = [
            {
                "label": "{} ({})".format(run["name"],
                                          len(run["samples"])),
                "value": run["name"]
            } for run in run_list]
        return [
            html.Div(
                [
                    html.Div(
                        dcc.Dropdown(
                            id="run-list",
                            options=run_list_options,
                            placeholder="Sequencing run"
                        )
                    )
                ],
                className="nine columns"
            ),
            html.Div(
                [
                    dcc.Link(
                        "Go to run",
                        id="run-link",
                        href="/",
                        className="button button-primary")
                ],
                className="three columns"
            )
        ]
    else:
        return


@app.callback(
    Output("run-link", "href"),
    [Input("run-list", "value")]
)
def update_run_button(run):
    if run is not None:
        return "/" + run
    else:
        return "/"


application = app.server  # Required for uwsgi

app.run_server(debug=True, host="0.0.0.0")
