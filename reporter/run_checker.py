# -*- coding: utf-8 -*-
import os
import sys

import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import plotly.figure_factory as ff
from dash.dependencies import Input, Output, State

import dash_scroll_up

import components.mongo_interface
import import_data


#Globals
PAGESIZE = 100
COMPONENTS = ['whats_my_species', 'qcquickie', 'assemblatron', 'analyzer', 'testomatic']

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
        html.Div(id="run-selector", className="row"),
        html.Div(id="run-report")

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

@app.callback(
    Output("run-report", "children"),
    [Input("run-name", "children")]
)
def update_run_report(run):
    data = []
    figure = {}
    if run != "":
        y = []
        x = COMPONENTS
        z = []
        z_text = []
        samples = import_data.get_sample_component_status(run)
        for name, s_components in samples.items():
            sample_z = []
            sample_z_text = []
            y.append(name)
            for component in COMPONENTS:
                if component in s_components.keys():
                    sample_z.append(s_components[component][0])
                    sample_z_text.append(s_components[component][1])
                else:
                    sample_z.append(0)
                    sample_z_text.append("None")
            z.append(sample_z)
            z_text.append(sample_z_text)
        # transpose z
        #z = [list(x) for x in zip(*z)]
        #z_text = [list(x) for x in zip(*z_text)]
        # print(x)
        # print(y)
        # print(z)
        figure = ff.create_annotated_heatmap(z=z,
                                             x=x,
                                             y=y,
                                             zmin=-1,
                                             zmax=2,
                                             annotation_text=z_text,
                                             colorscale=[[0, 'red'], [0.33, 'white'], [1, 'green']])
    
    figure.layout.margin.update({"l":200})
    # figure["layout"] = go.Layout(
    #     hovermode="closest",
    #     title="Run plot",
    #     margin=go.layout.Margin(
    #         l=175,
    #         r=50,
    #         b=125,
    #         t=50
    #     ),
    #     yaxis={"tickfont": {"size": 10}},
    #     xaxis={"showgrid": True}
    # )

    return [dcc.Graph(
        id="run-graph",
        figure=figure,
        style={"height": "1600px"}
    )]
    



application = app.server  # Required for uwsgi

#app.run_server(debug=True, host="0.0.0.0")
