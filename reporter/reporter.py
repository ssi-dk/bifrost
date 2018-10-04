# -*- coding: utf-8 -*-
import os
import sys
import urllib

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dt
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State

import dash_scroll_up

import flask  #used for the image server


import components.mongo_interface
import import_data
from components.table import html_table, html_td_percentage
from components.summary import html_div_summary
from components.sample_report import children_sample_list_report, generate_sample_folder
from components.images import list_of_images, static_image_route, get_species_color, COLOR_DICT, image_directory
import components.global_vars as global_vars

#Globals
#also defined in mongo_interface.py
PAGESIZE = 25


def hex_to_rgb(value):
    value = value.lstrip("#")
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def short_species(species):
    words = species.split(" ")
    if len(words) == 1:
        return species
    return "{}. {}".format(words[0][0], " ".join(words[1:]))




app = dash.Dash()
app.title = "Serum QC"
app.config["suppress_callback_exceptions"] = True

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
        html.Div(dt.DataTable(rows=[{}], editable=False), style={"display": "none"}),
        html.H1("SerumQC REPORT"),
        html.H2("", id="run-name"),
        html.H3(html.A("Wiki is here", href="https://teams.microsoft.com/l/channel/19%3a7b0b9a088602419e9f84630bacc84c2e%40thread.skype/tab%3a%3a9098abb1-75f5-410a-9011-87db7d42f3c2?label=Wiki&groupId=16852743-838a-400e-921d-6c50cc495b2f&tenantId=d0155445-8a4c-4780-9c13-33c78f22890e")),
        dcc.Location(id="url", refresh=False),
        html.Div(html_table([["run_name", ""]]), id="run-table"),
        html_div_summary(),
        html.Div(
            [
                html.H5(
                    [
                        "Update report (",
                        html.Span(id="report-count"),
                        " samples selected)"
                    ],
                    className="box-title"
                    ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Button(
                                    "QCQuickie",
                                    id="update-qcquickie",
                                    n_clicks_timestamp=0,
                                    className="button-primary u-full-width"
                                )
                            ],
                            className="three columns"
                        ),
                        html.Div(
                            [
                                html.Button(
                                    "Assemblatron",
                                    id="update-assemblatron",
                                    n_clicks_timestamp=0,
                                    className="button-primary u-full-width"
                                )
                            ],
                            className="three columns"
                        ),
                        html.Div(
                            [
                                html.Button(
                                    "Analyzer",
                                    id="update-analyzer",
                                    n_clicks_timestamp=0,
                                    className="button-primary u-full-width"
                                )
                            ],
                            className="three columns"
                        ),
                        html.Div(
                            [
                                html.Button(
                                    "Table report",
                                    id="update-table",
                                    n_clicks_timestamp=0,
                                    className="button-primary u-full-width"
                                )
                            ],
                            className="three columns"
                        )
                    ],
                    className="row",
                    style={"marginBottom": "15px"}
                ),
                html.Div([
                        html.Div(
                            [
                                html.Button(
                                    "Sample folder",
                                    id="generate-folder",
                                    n_clicks_timestamp=0,
                                    className="button-primary u-full-width"
                                )
                            ],
                            className="three columns"
                        )
                    ],
                    className="row"
                )
            ],
            className="border-box"
        ),
        html.Div(id="current-report"),
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
    Output("group-form", "className"),
    [Input("url", "pathname")]
)
def hide_group_if_in_url(pathname):
    if pathname is None:
        pathname = "/"
    path = pathname.split("/")
    if len(path) > 2 and path[2] != "":
        return "hidden"
    else:
        return ""

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
    Output("report-count", "children"),
    [Input("summary-plot", "selectedData"),
        Input("selected-samples-list", "value")]
)
def display_selected_data(plot_selected, selected_samples_list):
    if plot_selected is not None and len(plot_selected["points"]):
        return len([sample["text"]
                    for sample in plot_selected["points"]])
    else:
        return len(selected_samples_list[0].split("\n"))

@app.callback(
    Output("lasso-div", "children"),
    [Input("summary-plot", "selectedData"),
        Input("report-count", "children")]
)
def display_selected_data(selected_data, ignore_this):
    # ignore_this is there so the function is called 
    # when the sample list is updated.
    if selected_data is not None and len(selected_data["points"]):
        points = [sample["text"]
                    for sample in selected_data["points"]]
        sample_ids = [sample["customdata"]
                        for sample in selected_data["points"]] 
        return [
            html.Label(
                [
                    "Selected from plot lasso (",
                        str(len(points)),
                    "):"
                ],
                htmlFor="selected-from-plot"),
            dcc.Textarea(
                id="selected-from-plot",
                className="u-full-width",
                style={"resize": "none"},
                readOnly=True,
                value=", ".join(points)
            ),
            html.Div(",".join(sample_ids),
                        style={"display": "none"},
                        id="lasso-sample-ids")
        ]
    else:
        return [html.Div("",
                        style={"display": "none"},
                        id="lasso-sample-ids")]

@app.callback(
    Output("group-div", "children"),
    [Input("run-name", "children"),
    Input("url", "pathname")]
)
def update_group_list(run_name, pathname):
    if len(run_name) == 0:
        group_list = import_data.get_group_list()
    else:
        group_list = import_data.get_group_list(run_name)
    group_options = []
    group_list_options = []
    for item in group_list:
        if item["_id"] == None:
            group_options.append("Not defined")
            group_list_options.append({
                "label": "Not defined ({})".format(item["count"]),
                "value": "Not defined"
            })
        else:
            group_options.append(item["_id"])
            group_list_options.append({
                "label": "{} ({})".format(item["_id"], item["count"]),
                "value": item["_id"]
            })
    if pathname is None:
        pathname = "/"
    path = pathname.split("/")
    if len(path) > 2 and path[2] != "":
        group_options = [path[2]]
    return dcc.Dropdown(
        id="group-list",
        options=group_list_options,
        multi=True,
        value=group_options
    )

@app.callback(
    Output("species-div", "children"),
    [Input("run-name", "children")]
)
def update_species_list(run_name):
    if len(run_name) == 0:
        species_list = import_data.get_species_list()
    else:
        species_list = import_data.get_species_list(run_name)
    species_options = []
    species_list_options = []
    for item in species_list:
        if item["_id"] == None:
            species_options.append("Not classified")
            species_list_options.append({
                "label": "Not classified ({})".format(item["count"]),
                "value": "Not classified"
            })
        else:
            species_options.append(item["_id"])
            species_list_options.append({
                "label": "{} ({})".format(short_species(item["_id"]), item["count"]),
                "value": item["_id"]
            })
    return dcc.Dropdown(
        id="species-list",
        options=species_list_options,
        multi=True,
        value=species_options
    )


@app.callback(
    Output("run-table", "children"),
    [Input("run-name", "children"),
        Input("url", "pathname")]
)
def update_run_table(run_name, pathname):
    if pathname is None:
        pathname = "/"
    path = pathname.split("/")
    if run_name == "Not found":
        run = "Run not found!"
    elif run_name == None or run_name == "":
        run = "No run selected"
    else:
        run = run_name
    if len(path) > 2 and path[2] != "":
        group = path[2]
        return html_table([["Run Name", run], ["Supplying lab", group]])
    else:
        return html_table([["Run Name", run]])

@app.callback(
    Output(component_id="page-n",
            component_property="children"),
    [Input(component_id="prevpage", component_property="n_clicks_timestamp"),
        Input(component_id="nextpage", component_property="n_clicks_timestamp")],
    [State(component_id="page-n", component_property="children"),
        State(component_id="max-page", component_property="children")]
)
def next_page(prev_ts, next_ts, page_n, max_page):
    page_n = int(page_n)
    max_page = int(max_page)
    if prev_ts > next_ts:
        return max(page_n - 1, 0)
    elif next_ts > prev_ts:
        return min(page_n + 1, max_page)
    else:
        return 0


@app.callback(
    Output(component_id="sample-report", component_property="children"),
    [Input(component_id="page-n", component_property="children"),
        Input(component_id="sample-report", component_property="data-content")],
    [State("lasso-sample-ids", "children"),
        State("selected-samples-ids", "children")]
        )
def sample_report(page_n, data_content, lasso_selected, prefilter_samples):
    if data_content not in ["qcquickie", "assemblatron", "analyzer"]:
        return []
    page_n = int(page_n)
    if lasso_selected != "":
        samples = lasso_selected.split(",")  # lasso first
    else:
        samples = prefilter_samples.split(",")
    #NOTE Could optimize this by not getting all sample's info from mongo before paginating
    page = import_data.filter_all(sample_ids=samples, page=page_n)
    page = page.sort_values(["species","name"])
    page = import_data.add_sample_runs(page)
    max_page = len(samples) // PAGESIZE
    page_species = page["species"].unique().tolist()
    species_plot_data = import_data.get_species_plot_data(page_species, page["_id"].tolist())
    return [
        html.H4("Page {} of {}".format(page_n + 1, max_page + 1)),
        html.Div(children_sample_list_report(page, data_content, species_plot_data))
    ]

@app.callback(
    Output(component_id="prevpage", component_property="disabled"),
    [Input(component_id="page-n", component_property="children")]
)
def update_prevpage(page_n):
    if int(page_n) == 0:
        return True
    else:
        return False

@app.callback(
    Output(component_id="nextpage", component_property="disabled"),
    [Input(component_id="page-n", component_property="children"),
        Input(component_id="max-page", component_property="children")]
)
def update_nextpage(page_n, max_page):
    page_n = int(page_n)
    max_page = int(max_page)
    if page_n == max_page:
        return True
    else:
        return False

@app.callback(
    Output("current-report", "children"),
    [Input("update-qcquickie", "n_clicks_timestamp"),
        Input("update-assemblatron", "n_clicks_timestamp"),
        Input("update-analyzer", "n_clicks_timestamp"),
        Input("update-table", "n_clicks_timestamp"),
        Input("generate-folder", "n_clicks_timestamp")],
    [State("lasso-sample-ids", "children"),
        State("selected-samples-ids", "children")]
)
def update_report(n_qcquickie_ts, n_assemblatron_ts,
                n_analyzer_ts, n_table_ts, n_generate_ts,
                lasso_selected, prefilter_samples):
    if lasso_selected != "":
        samples = lasso_selected.split(",")  # lasso first
    elif prefilter_samples != "":
        samples = prefilter_samples.split(",")
    else:
        samples = []
    
    max_page = len(samples) // PAGESIZE

    last_module_ts = max(
        n_qcquickie_ts, n_assemblatron_ts, n_analyzer_ts, n_generate_ts, n_table_ts)
    if min(n_qcquickie_ts, n_assemblatron_ts,
            n_analyzer_ts, n_generate_ts, n_table_ts) == last_module_ts:
        return []
    report = False
    if n_qcquickie_ts == last_module_ts:
        title = "QCQuickie Report"
        content = "qcquickie"
        report = True
    elif n_assemblatron_ts == last_module_ts:
        title = "Assemblatron Report"
        content = "assemblatron"
        report = True
    elif n_analyzer_ts == last_module_ts:
        title = "Analyzer Report"
        content = "analyzer"
        report = True
    if report:
        return [
            html.H3(title),
            html.Span("0", style={"display": "none"}, id="page-n"),
            html.Span(max_page, style={"display": "none"}, id="max-page"),
            html.Div(
                [
                    html.Div(
                        [
                            html.Button(
                                "Previous page", id="prevpage", n_clicks_timestamp=0),
                        ],
                        className="three columns"
                    ),
                    html.Div(
                        [
                            # html.H4("Page {} of {}".format(page_n + 1, max_page + 1))
                        ],
                        className="three columns"
                    ),
                    html.Div(
                        [
                            html.Button(
                                "Next page", id="nextpage", n_clicks_timestamp=0),
                        ],
                        className="three columns"
                    ),
                ],
                className="row"
            ),
            
            html.Div(id="sample-report", **{"data-content": content}),
        ]
    elif n_table_ts == last_module_ts:  # table was clicked
        dataframe = import_data.filter_all(sample_ids=samples)
        if "analyzer.ariba_resfinder" in dataframe:
            dataframe = dataframe.drop(
            columns="analyzer.ariba_resfinder")
        csv_string = dataframe.to_csv(index=False, encoding="utf-8", sep="\t")
        csv_string = 'data:text/tab-separated-values;charset=utf-8,' + urllib.parse.quote(csv_string)
        return [
            html.H3("Table Report"),
            # We have to drop those columns with nested values for the table
            html.A("Download run (tsv)", href=csv_string, download='report.tsv'),
            dt.DataTable(
                rows=dataframe.to_dict("records"),

                # columns=global_vars.columns, # sets the order
                # column_widths=[150]*len(global_vars.columns),
                editable=False,
                filterable=True,
                sortable=True,
                selected_row_indices=[],
                id="datatable-samples"
            )
        ]
    elif n_generate_ts == last_module_ts:
        return generate_sample_folder(samples)
    return []

@app.callback(
    Output(component_id="selected-samples", component_property="children"),
    [Input(component_id="species-list", component_property="value"),
        Input(component_id="group-list", component_property="value"),
        Input(component_id="run-name", component_property="children")]
)
def update_selected_samples(species_list, group_list, run_name):

    samples = import_data.filter_name(species=species_list,
        group=group_list, run_name=run_name)
    sample_names = []
    sample_ids = []
    for sample in samples:
        try:
            sample_names.append("{}".format(sample["name"]))
        except KeyError:
            sample_names.append(sample["sample_sheet"]["sample_name"])
        sample_ids.append(str(sample["_id"]))
    return [
        html.Label(
            [
                "Selected Samples (",
                len(samples),
                "):"
            ],
            htmlFor="plot-list"),
        dcc.Textarea(
            className="u-full-width",
            style={"resize": "none",
                    "height": "300px",
                    "white-space": "pre"},
            readOnly=True,
            value=["\n".join(sample_names)],
            id="selected-samples-list"
        ),
        html.Div(",".join(sample_ids),
                    style={"display": "none"},
                    id="selected-samples-ids")
    ]


@app.callback(
    Output(component_id="testomatic-report", component_property="children"),
    [Input(component_id="species-list", component_property="value"),
        Input(component_id="group-list", component_property="value"),
        Input(component_id="run-name", component_property="children")]
)
def update_test_table(species_list, group_list, run_name):

    columns = ["name", "species", "sample_sheet.group", "sample_sheet.Comments", 'testomatic.qcquickie.action',
            'testomatic.assemblatron.action', 'testomatic.assemblatron.10xgenomesize',
            'testomatic.assemblatron.1x25xsizediff',
            'testomatic.assemblatron.1xgenomesize', 'testomatic.assemblatron.avgcoverage',
            'testomatic.assemblatron.numreads', 'testomatic.base.readspresent',
            'testomatic.qcquickie.10xgenomesize',
            'testomatic.qcquickie.1x25xsizediff',
            'testomatic.qcquickie.1xgenomesize',
            'testomatic.qcquickie.avgcoverage', 'testomatic.qcquickie.numreads',
            'testomatic.whats_my_species.maxunclassified',
            'testomatic.whats_my_species.minspecies',
            'testomatic.whats_my_species.nosubmitted',
            'testomatic.whats_my_species.submitted==detected', "_id"]
    column_names = ['name', 'Species', 'supplying lab', "Comments", 'qcquickie QC',
                    'assemblatron QC', 'assemblatron 10xgenomesize',
                    'assemblatron 1x25xsizediff', 'assemblatron 1xgenomesize',
                    'assemblatron avgcoverage', 'assemblatron numreads', 'base readspresent',
                    'qcquickie 10xgenomesize', 'qcquickie 1x25xsizediff',
                    'qcquickie 1xgenomesize', 'qcquickie avgcoverage',
                    'qcquickie numreads', 'whats_my_species maxunclassified',
                    'whats_my_species minspecies', 'whats_my_species nosubmitted',
                    'whats_my_species submitted==detected', '_id']
    samples_df = import_data.filter_all(
        species_list, group_list, run_name)

    samples_df.species = list(map(lambda x:short_species(x), samples_df.species))
    table = [{"list":column_names,"className":"header"}]
    for sample in samples_df.iterrows():
        row = []
        for value, column in zip(sample[1][columns], columns):
            if value.startswith("fail") or value.startswith("undefined") \
            or value == "supplying lab" or value =="core facility":
                td = html.Td(str(value), className="cell red")
            else:
                td = html.Td(str(value), className="cell")
            row.append(td)
        table.append(html.Tr(row, className="trow"))

    return [html.Table(table)]


@app.callback(
    Output(component_id="summary-plot", component_property="figure"),
    [Input(component_id="species-list", component_property="value"),
        Input(component_id="group-list", component_property="value"),
        Input(component_id="run-name", component_property="children"),
        Input(component_id="plot-list", component_property="value")]
)
def update_coverage_figure(species_list, group_list, run_name, plot_value):
    plot_query = global_vars.PLOTS[plot_value]["projection"]
    plot_func = global_vars.PLOTS[plot_value].get("func")
    
    data = []
    plot_df = import_data.filter_all(species_list, group_list, run_name, plot_func)
    if species_list is None: species_list = []
    species_count = 0
    if 'species' in plot_df:
        # reverse the list so it looks right on plot
        for species in reversed(species_list):
            species_df = plot_df[plot_df.species == species]
            if species == "Not classified":
                species_name = species
            else:
                species_name = "<i>{}</i>".format(short_species(species))
            if (plot_query in species_df):
                if (len(species_df)): species_count += 1
                data.append(
                    go.Box(
                        x=species_df.loc[:, plot_query],
                        text=species_df["name"],
                        marker=dict(
                            color=COLOR_DICT.get(species, None),
                            size=4
                        ),
                        
                        boxpoints="all",
                        jitter=0.3,
                        pointpos=-1.8,
                        name="{} ({})".format(species_name,species_df["_id"].count()),
                        showlegend=False,
                        customdata=species_df["_id"]
                    )
            )
    height = max(450, species_count*20 + 200)
    return {
        "data": data,
        "layout": go.Layout(
            hovermode="closest",
            title=plot_value.replace("_", " "),
            height=height,
            margin=go.layout.Margin(
                l=175,
                r=50,
                b=25,
                t=50
            ),
            yaxis={"tickfont":{"size": 10}},
            xaxis={"showgrid": True}
        )
    }

@app.callback(
    Output("group-list", "value"),
    [Input("group-all", "n_clicks"),
        Input("url", "pathname")],
    [State("run-name", "children")]
)
def all_groups(n_clicks, pathname, run_name):
    if pathname is None:
        pathname = "/"
    path = pathname.split("/")
    if len(path) > 2 and path[2] != "":
        return [path[2]]
    if len(run_name) == 0:
        group_list = import_data.get_group_list()
    else:
        group_list = import_data.get_group_list(run_name)
    group_options = []
    for item in group_list:
        if item["_id"] == None:
            group_options.append("Not defined")
        else:
            group_options.append(item["_id"])

    return group_options

@app.callback(
    Output("species-list", "value"),
    [Input("species-all", "n_clicks")],
    [State("run-name", "children")]
)
def all_species(n_clicks, run_name):
    if len(run_name) == 0:
        species_list = import_data.get_species_list()
    else:
        species_list = import_data.get_species_list(run_name)
    species_options = []
    for item in species_list:
        if item["_id"] == None:
            species_options.append("Not classified")
        else:
            species_options.append(item["_id"])

    return species_options

application = app.server # Required for uwsgi

