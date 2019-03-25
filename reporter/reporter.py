# -*- coding: utf-8 -*-
import os
import sys
import urllib
import datetime
import time
from io import StringIO

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import dash_auth
import pandas as pd
import numpy as np
import plotly.graph_objs as go
from plotly import tools 
from dash.dependencies import Input, Output, State

from flask import request   # To get client IP for pass/fail stamp

import dash_scroll_up

import keys

import components.mongo_interface
import components.import_data as import_data
from components.table import html_table, html_td_percentage
from components.summary import html_div_summary, filter_notice_div
from components.sample_report import children_sample_list_report, generate_sample_folder
from components.images import list_of_images, static_image_route, image_directory
import components.global_vars as global_vars
import components.admin as admin

# Globals
# also defined in mongo_interface.py
PAGESIZE = 25
ADMIN = False


def hex_to_rgb(value):
    value = value.lstrip("#")
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def short_species(species):
    if species is None or pd.isna(species):
        return None
    words = species.split(" ")
    if len(words) == 1:
        return species
    return "{}. {}".format(words[0][0], " ".join(words[1:]))


app = dash.Dash()
app.title = "bifrost"
app.config["suppress_callback_exceptions"] = True

if hasattr(keys, "pass_protected") and keys.pass_protected:
    dash_auth.BasicAuth(
        app,
        keys.USERNAME_PASSWORD
    )

# Temp css to make it look nice
# Dash CSS
app.css.append_css(
    {"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})
# Lato font
app.css.append_css(
    {"external_url": "https://fonts.googleapis.com/css?family=Lato"})

app.layout = html.Div([
    html.Div(className="container", children=[
        html.Div(id="placeholder0", style={"display": "none"}),
        html.Div(id="placeholder1", style={"display": "none"}),
        html.Div(id="placeholder2", style={"display": "none"}),
        dash_scroll_up.DashScrollUp(
            id="input",
            label="UP",
            className="button button-primary no-print"
        ),
        html.Div(dash_table.DataTable(editable=False), style={"display": "none"}),
        html.H1("bifrost report", style={"display": "none"}),
        html.Img(
            src="/assets/img/report.png",
            className="main-logo"
        ),
        html.H2("Loading...", id="run-name", style={"display": "none"}),
        html.Div([
            html.Div(id="report-link", className="u-pull-left"),
            html.H6(html.A("Wiki", href="https://teams.microsoft.com/l/channel/19%3a7b0b9a088602419e9f84630bacc84c2e%40thread.skype/tab%3a%3a9098abb1-75f5-410a-9011-87db7d42f3c2?label=Wiki&groupId=16852743-838a-400e-921d-6c50cc495b2f&tenantId=d0155445-8a4c-4780-9c13-33c78f22890e"), className="u-pull-right"),
        ], className="row"),
        html.Details([
            html.Summary("Latest changes..."),
            html.Div([
                html.P(
                    "2019-03-19: Add Priority to run checker."),
                html.P(
                    "2019-03-14: Add Supplying lab feedback option."),
                html.P(
                    "2019-03-07: Add QC values to Run checker."),
                html.P(
                    "2019-02-28: Add assemblatron links to sample folder generator. Only samples run from now on will have the right paths."),
                html.P("2019-02-13: Filters default to empty (but behavior is the same), bring back apply button to reduce performance issues, add current run to run selector and other performance fixes. Add Latest changes section."),
            ])
        ]),
        dcc.Store(id="data-store", storage_type="memory"),
        dcc.Location(id="url", refresh=False),
        html.Div(html_table([["run_name", ""]]), id="run-table"),
        html_div_summary(),
        html.Div(id="current-report"),
        html.Div(
            [
                html.H5(
                    [
                        "Generate sample folder script"
                    ],
                    className="box-title"
                    ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Button(
                                    "Assembly/Sample folder",
                                    id="generate-folder",
                                    n_clicks_timestamp=0,
                                    className="button-primary u-full-width"
                                )
                            ],
                            className="twelve columns"
                        ),
                    ],
                    className="row"
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(id="sample-folder-div")
                            ],
                            className="twelve columns"
                        ),
                    ],
                    className="row",
                    style={"marginBottom": "15px"}
                )
            ],
            className="border-box"
        ),
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
    Output("group-list", "value"),
    [Input("url", "pathname")]
)
def hide_group_if_in_url(pathname):
    if pathname is None:
        pathname = "/"
    path = pathname.split("/")
    if len(path) > 2 and path[2] != "":
        return [path[2]]
    else:
        return []


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
    Output("lasso-div", "children"),
    [Input("summary-plot", "selectedData"),
     Input('datatable-ssi_stamper', 'derived_virtual_data'),
     Input('datatable-ssi_stamper', 'derived_virtual_selected_rows')]
)
def display_selected_data(selected_data, rows, selected_rows):
    # ignore_this is there so the function is called 
    # when the sample list is updated.
    if rows == [{}] or rows == [] or rows == None:
        return [
            html.Label(
                [
                    "Selected from table/lasso (0):"
                ],
                htmlFor="selected-from-plot"),
            dcc.Textarea(
                id="selected-from-plot",
                className="u-full-width",
                style={"resize": "none"},
                readOnly=True,
                value=""
            ),
            html.Div(style={"display": "none"},
                        id="lasso-sample-ids")
        ]
    dtdf = pd.DataFrame(rows)
    if selected_rows is not None and len(selected_rows) > 0:
        dtdf = dtdf.iloc[selected_rows]
    points = list(map(str, list(dtdf["name"])))
    sample_ids = list(dtdf["_id"])
    if selected_data is not None and len(selected_data["points"]):
        lasso_points = set([sample["text"]
                    for sample in selected_data["points"]])
        lasso_sample_ids = set([sample["customdata"]
                        for sample in selected_data["points"]])
        union_points = set(points).intersection(lasso_points)
        union_sample_ids = set(sample_ids).intersection(lasso_sample_ids)
        # This way we keep the table order.
        points = list(dtdf[dtdf["_id"].isin(lasso_sample_ids)]["name"])
        sample_ids = [sample_id for sample_id in sample_ids if sample_id in union_sample_ids]
    return [
        html.Label(
            [
                "Selected from table/lasso (",
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


@app.callback(
    Output("group-div", "children"),
    [Input("run-name", "children")],
    [State("url", "pathname")]
)
def update_group_list(run_name, pathname):
    if run_name == "Loading...":
        return None
    if len(run_name) == 0:
        group_list = import_data.get_group_list()
    else:
        group_list = import_data.get_group_list(run_name)
    group_list_options = []
    for item in group_list:
        if item["_id"] == None:
            group_list_options.append({
                "label": "Not defined ({})".format(item["count"]),
                "value": "Not defined"
            })
        else:
            group_list_options.append({
                "label": "{} ({})".format(item["_id"], item["count"]),
                "value": item["_id"]
            })
    if pathname is None:
        pathname = "/"
    path = pathname.split("/")
    if len(path) > 2 and path[2] != "":
        value = [path[2]]
    else:
        value = []
    return dcc.Dropdown(
        id="group-list",
        options=group_list_options,
        multi=True,
        placeholder="All groups selected",
        value=value
    )

@app.callback(
    Output("species-div", "children"),
    [Input("run-name", "children"),
    Input("form-species-source", "value")]
)
def update_species_list(run_name, form_species):
    if run_name == "Loading...":
        return None
    if len(run_name) == 0:
        species_list = import_data.get_species_list(form_species)
    else:
        species_list = import_data.get_species_list(form_species, run_name)
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
        placeholder="All species selected",
        multi=True,
        value=[]
    )


@app.callback(
    Output("run-list-div", "children"),
    [Input("run-name", "children")]
)
def update_run_options(run_name):
    run_list = import_data.get_run_list()
    options = [
        {
            "label": "{} ({})".format(run["name"],
                                      len(run["samples"])),
            "value": run["name"]
        } for run in run_list]
    if run_name not in ["", "Loading...", "Not Found"]:
        return dcc.Dropdown(
            id="run-list",
            options=options,
            value=run_name
        )
    else:
        return dcc.Dropdown(
            id="run-list",
            options=options,
            placeholder="Sequencing run"
        )


@app.callback(
    Output("run-table", "children"),
    [Input("run-name", "children"),
        Input("url", "pathname")]
)
def update_run_table(run_name, pathname):
    if run_name == "Loading...":
        return None
    if pathname is None:
        pathname = "/"
    path = pathname.split("/")
    if run_name == "Not found":
        run = "Run not found!"
    elif run_name == None or run_name == "":
        run = "No run selected"
    else:
        run = run_name
        year = run[:2]
        run_path = keys.run_path + "20{}/{}".format(year, run_name)
        run_link = "file:/" + run_path
        if hasattr(keys, "path_platform_windows") and keys.path_platform_windows:
            run_path = run_path.replace("/", "\\")
        run_a=html.A(run_path, href=run_link)
        run_checker_url = "{}/{}".format(keys.run_checker_url, run_name)
        run_checker_link = html.A(run_checker_url, href=run_checker_url)
        if len(path) > 2 and path[2] != "":
            group = path[2]
            return html_table([
                ["Run Name", run],
                ["Run Path", run_a],
                ["Supplying lab", group],
                ["Run Checker", run_checker_link]
            ])
        else:
            return html_table([
                ["Run Name", run],
                ["Run Path", run_a],
                ["Run Checker", run_checker_link]
            ])
    return html_table([["Run Name", run]])

@app.callback(
    Output("page-n",
            "children"),
    [Input("prevpage", "n_clicks_timestamp"),
        Input("prevpage2", "n_clicks_timestamp"),
        Input("nextpage", "n_clicks_timestamp"),
        Input("nextpage2", "n_clicks_timestamp")],
    [State("page-n", "children"),
        State("max-page", "children")]
)
def next_page(prev_ts, prev_ts2, next_ts, next_ts2, page_n, max_page):
    page_n = int(page_n)
    max_page = int(max_page)
    if max(prev_ts, prev_ts2) > max(next_ts, next_ts2):
        return max(page_n - 1, 0)
    elif max(next_ts, next_ts2) > max(prev_ts, prev_ts2):
        return min(page_n + 1, max_page)
    else:
        return 0


@app.callback(
    Output("sample-report", "children"),
    [Input("page-n", "children")],
    [State("lasso-sample-ids", "children"),
        State("data-store", "data")]
        )
def sample_report(page_n, lasso_selected, data_store):
    page_n = int(page_n)
    #json_data = StringIO(data_store)
    data = pd.DataFrame.from_dict(data_store)
    if lasso_selected != "" and lasso_selected is not None:
        lasso = lasso_selected.split(",")  # lasso first
        data = data[data._id.isin(lasso)]
    if len(data) == 0: return []
    samples = data["_id"]
    data = data.sort_values(["species","name"])
    skips = PAGESIZE * (page_n)
    page = data[skips:skips+PAGESIZE]
    page = import_data.add_sample_runs(page)
    max_page = len(samples) // PAGESIZE
    page_species = page["species"].unique().tolist()
    species_plot_data = import_data.get_species_plot_data(page_species, page["_id"].tolist())
    # We need to have fake radio buttons with the same ids to account for times 
    # when not all PAGESIZE samples are shown and are not taking the ids required by the callback
    html_fake_radio_buttons = html.Div([dcc.RadioItems(
        options=[
            {'label': '', 'value': 'nosample'}
        ],
        value='noaction',
        id="sample-radio-{}".format(n_sample)
    ) for n_sample in range(len(page), PAGESIZE)], style={"display": "none"})
    return [
        html.H4("Page {} of {}".format(page_n + 1, max_page + 1)),
        html.Div(children_sample_list_report(page, species_plot_data)),
        html_fake_radio_buttons,
        admin.html_qc_expert_form(),
        html.H4("Page {} of {}".format(page_n + 1, max_page + 1)),
        html.Div(id="placeholder0", style={"display": "none"}),
        dcc.ConfirmDialog(
            id='qc-confirm',
            message='Are you sure you want to send sample feedback?',
        )
    ]

@app.callback(
    Output("sample-folder-div", "children"),
    [
        Input("generate-folder", "n_clicks_timestamp")],
    [
        State("lasso-sample-ids", "children"),
        State("data-store", "data")]
)
def generate_sample_folder_div(n_generate_ts,
                  lasso_selected, data_store):

    if n_generate_ts == 0:
        return None

    if lasso_selected != "":
        samples = lasso_selected.split(",")  # lasso first
    elif data_store != None:
        #json_data = StringIO(data_store)
        samples = pd.DataFrame.from_dict(data_store)["_id"]
    else:
        samples = []

    return generate_sample_folder(samples)


@app.callback(
    Output("current-report", "children"),
    [Input("lasso-sample-ids", "children")]
)

def update_report(lasso_selected):
    if lasso_selected is None or lasso_selected == "":
        return []
    else:
        sample_n = lasso_selected.count(",") + 1
        return [
            html.H3("Sample Report"),
            html.Span("0", style={"display": "none"}, id="page-n"),
            html.Span(sample_n // PAGESIZE,
                    style={"display": "none"}, id="max-page"),
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
            
            html.Div(id="sample-report"),

            html.Div(
                [
                    html.Div(
                        [
                            html.Button(
                                "Previous page", id="prevpage2", n_clicks_timestamp=0),
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
                                "Next page", id="nextpage2", n_clicks_timestamp=0),
                        ],
                        className="three columns"
                    ),
                ],
                className="row"
            ),
        ]


@app.callback(
    Output("apply-filter-button", "children"),
    [Input("species-list", "value"),
        Input("group-list", "value"),
        Input("qc-list", "value"),
        Input("run-name", "children"),
        Input("apply-filter-button", "n_clicks_timestamp")]
)
def update_filter_ts(species_list, group_list, qc_list, run_name, button_timestamp):
    button_ts = button_timestamp / 1000
    server_ts = datetime.datetime.now().timestamp()
    if button_ts == 0 or abs(button_ts - server_ts) < 1:
        return "Apply filter"
    return "Apply filter (click to reload)"

@app.callback(
    Output("data-store", "data"),
    [Input("apply-filter-button", "n_clicks_timestamp"),
     Input("run-name", "children")],
    [State("species-list", "value"),
     State("form-species-source", "value"),
        State("group-list", "value"),
        State("qc-list", "value")]
)
def update_selected_samples(apply_button_ts, run_name, species_list, species_source, group_list, qc_list):
    if run_name == "Loading..." or \
        None in (species_list, group_list, qc_list, run_name):
        return '""'
    else:
        samples = import_data.filter_all(species=species_list, species_source=species_source,
            group=group_list, qc_list=qc_list, run_name=run_name)
    return samples.to_dict()


@app.callback(
    Output("ssi_stamper-report",
           "children"), 
    [Input("data-store", "data")]
)
def update_test_table(data_store):
    empty_table = [
        html.H6('No samples loaded. Click "Apply Filter" to load samples.'),
        filter_notice_div,
        html.Div([
            html.Div([
                "Download Table ",
                html.A("(tsv, US format)",
                       download='report.tsv'),
                " - ",
                html.A("(csv, EUR Excel format)",
                       download='report.csv')
            ], className="six columns"),
        ], className="row"),
        html.Div(dash_table.DataTable(id="datatable-ssi_stamper",
                                      data=[{}]), style={"display": "none"})
    ]
    if data_store == '""':
        return empty_table
    ##json_data = StringIO(data_store)
    tests_df = pd.DataFrame.from_dict(data_store)
    if len(tests_df) == 0:
        return empty_table
    qc_action = "ssi_stamper.assemblatron:action"
    qc_action = "stamp.ssi_stamper.value"
    if qc_action not in tests_df:
        tests_df[qc_action] = np.nan
    else:
        tests_df[qc_action] = tests_df[qc_action].str.split(":", expand=True)[1]
 
    if "R1" not in tests_df:
        tests_df["R1"] = np.nan

    no_reads_mask = tests_df["R1"] == ""
    tests_df.loc[no_reads_mask, qc_action] = "core facility (no reads)"
    mask = pd.isnull(tests_df[qc_action])
    tests_df.loc[mask, qc_action] = "not tested"
    slmask = tests_df[qc_action] == "supplying lab"
    tests_df.loc[slmask, qc_action] = "warning: supplying lab"
    
    user_stamp_col = "stamp.supplying_lab_check.value"
    # Overload user stamp to ssi_stamper
    if user_stamp_col in tests_df.columns:
        user_OK_mask = tests_df[user_stamp_col] == "pass:OK"
        tests_df.loc[user_OK_mask, qc_action] = "*OK"
        user_sl_mask = tests_df[user_stamp_col] == "fail:supplying lab"
        tests_df.loc[user_sl_mask, qc_action] = "*warning: supplying lab"
        user_cf_mask = tests_df[user_stamp_col] == "fail:core facility"
        tests_df.loc[user_cf_mask, qc_action] = "*core facility"

    # Split test columns
    columns = tests_df.columns
    split_columns = [
        "ssi_stamper.assemblatron:1x10xsizediff",
        "ssi_stamper.whats_my_species:minspecies",
        "ssi_stamper.whats_my_species:nosubmitted",
        "ssi_stamper.whats_my_species:detectedspeciesmismatch"
    ]
    i = 0
    for column in columns:
        if column in split_columns:
            new = tests_df[column].str.split(":", expand=True)
            loc = tests_df.columns.get_loc(column)
            #tests_df.drop(columns = [column], inplace=True)
            tests_df.insert(loc, column + "_QC", new[0])
            tests_df.insert(loc + 1, column + "_text", new[2])
        i += 1

    if "analyzer.mlst_report" in columns:
        first_split = tests_df["analyzer.mlst_report"].str.split(
            ",", n=1, expand=True)
        if len(first_split.columns) == 2:
            second_split = first_split[0].str.split(":", n=1, expand=True)
            if len(second_split.columns) == 2:
                keyerrormask = second_split[1] == " 'ariba_mlst/mlst_report_tsv'"
                second_split.loc[keyerrormask, 1] = np.nan
                tests_df["analyzer_mlst_type"] = second_split[1]
                tests_df["analyzer_mlst_alleles"] = first_split[1]

    test_cols = [col for col in columns if (col.startswith(
        "ssi_stamper") and not col.startswith("ssi_stamper.qcquickie"))]
    def concatenate_failed(row):
        res = []
        for col in test_cols:
            test_name = col.split(":")[-1]
            if type(row[col]) == str:
                fields = row[col].split(":")
                if fields[0] in ["fail", "undefined"]:
                    res.append("Test {}: {}, {}".format(test_name, fields[0], fields[1]))
        row["ssi_stamper_failed_tests"] = ". ".join(res)
        return row
    
    # Round columns:
    for col in global_vars.ROUND_COLUMNS:
        if col in tests_df.columns:
            tests_df[col] = round(tests_df[col],3)

    tests_df = tests_df.apply(concatenate_failed, axis="columns")


    COLUMNS = global_vars.COLUMNS

    # HORRIBLE HACK: but must be done because of bug https://github.com/plotly/dash-table/issues/224
    # Delete as soon as filter works with column ids

    def simplify_name(name):
        return name.replace(":", "_").replace(".", "_").replace("=", "_")
        
    tests_df.columns = list(map(simplify_name, tests_df.columns))
    COLUMNS = []
    for column in global_vars.COLUMNS:
        column["id"] = simplify_name(column["id"])
        COLUMNS.append(column)

    # END OF HORRIBLE HACK

    # Generate conditional formatting:
    style_data_conditional = []
    conditional_columns = [ col for col in tests_df.columns if col.startswith("QC_") ]
    
    for status, color in ("fail", "#ea6153"), ("undefined", "#f1c40f"):
        style_data_conditional += list(map(lambda x: {"if": {
            "column_id": x, "filter": '{} eq "{}"'.format(x, status)}, "backgroundColor": color}, conditional_columns))
    
    for status, color in ("core facility", "#ea6153"), ("warning: supplying lab", "#f1c40f"):
        style_data_conditional += [{"if": {
            "column_id": qc_action, "filter": 'QC_action eq "{}"'.format(status)}, "backgroundColor": color}]


    table = dash_table.DataTable(

        data=tests_df.to_dict("rows"),
        style_table={
            'overflowX': 'scroll',
            'overflowY': 'scroll',
            'maxHeight': '480'
        },
        columns=COLUMNS,
        # n_fixed_columns=1,
        style_cell={
            'width': '200px',
            'padding': '0 15px'
        },
        style_cell_conditional=[
            {
                "if": {"column_id": "ssi_stamper_failed_tests"},
                "textAlign": "left"
            }
        ],
        n_fixed_rows=1,
        row_selectable="multi",
        filtering=True, #Front end filtering
        sorting=True,
        selected_rows=[],
        style_data_conditional=style_data_conditional,
        id="datatable-ssi_stamper"
    )

    rename_dict = {item["id"]:item["name"] for item in COLUMNS}

    renamed = tests_df.rename(rename_dict, axis='columns')#[

    missing_columns = [a for a in list(rename_dict.values()) if not a in list(renamed.columns)]

    # add missing columns
    for column in missing_columns:
        renamed[column] = np.nan
    
    # reorder columns
    renamed = renamed[list(rename_dict.values())]

    csv_string_eur = renamed.to_csv(index=False, encoding="utf-8", sep=";", decimal=",")
    tsv_string_us = renamed.to_csv(index=False, encoding="utf-8", sep="\t")
    full_csv_string_eur = 'data:text/csv;charset=utf-8,' + \
        urllib.parse.quote(csv_string_eur)
    full_tsv_string_us = 'data:text/tab-separated-values;charset=utf-8,' + \
        urllib.parse.quote(tsv_string_us)
    return [
        html.H6("Filtered samples ({}):".format(len(tests_df["_id"]))),
        filter_notice_div,
        html.Div([
            html.Div([
                "Download Table ",
                html.A("(tsv, US format)",
                       href=full_tsv_string_us,
                       download='report.tsv'),
                " - ",
                html.A("(csv, EUR Excel format)",
                       href=full_csv_string_eur,
                       download='report.csv')
            ], className="six columns"),
        ], className="row"),
        html.Div(table)
        ]


@app.callback(
    Output("plot-species-div", "children"),
    [Input('datatable-ssi_stamper', 'derived_virtual_data'),
     Input('datatable-ssi_stamper', 'derived_virtual_selected_rows'),
     Input("plot-species-source", "value")],
    [State("plot-species", "value")]
)
def plot_species_dropdown(rows, selected_rows, plot_species, selected_species):
    plot_df = pd.DataFrame(rows)
    if selected_rows is not None and len(selected_rows) > 0:
        plot_df = plot_df.iloc[selected_rows]
    # part of the HACK, replace with "properties.detected_species" when issue is solved
    species_col = "properties_detected_species"

    if plot_species == "provided":
        species_col = "properties_provided_species"
    elif plot_species == "detected":
        species_col = "properties_detected_species"
    # end HACK
    if species_col not in plot_df or plot_df[species_col].unique() is None:
        return dcc.Dropdown(
            id="plot-species"
        )
    plot_df.loc[pd.isnull(plot_df[species_col]), species_col] = "Not classified"
    species_list = plot_df[species_col].unique()
    species_list = ["All species",] + list(species_list)
    if selected_species == "Not classified" or selected_species is None or selected_species not in species_list:
        if species_list[0] == "Not classified" and len(species_list) > 1:
            selected_species = species_list[1]
        else:
            selected_species = species_list[0]
    species_list_options = [
        {
            "label": species,
            "value": species
        } for species in species_list]
    return dcc.Dropdown(
        id="plot-species",
        options=species_list_options,
        value=selected_species
    )


@app.callback(
    Output("summary-plot", "selectedData"),
    [Input("data-store", "data"),
     Input("plot-species", "value"),
     Input('datatable-ssi_stamper', 'derived_virtual_data'),
     Input('datatable-ssi_stamper', 'derived_virtual_selected_rows')]
)
def reset_selection(sample_ids, plot_value, rows, selected_rows):
    return {"points":[]}


@app.callback(
    Output("summary-plot", "figure"),
    [Input("plot-species", "value"),
     Input('datatable-ssi_stamper', 'derived_virtual_data'),
     Input('datatable-ssi_stamper', 'derived_virtual_selected_rows')],
    [State("plot-species-source", "value")]
)
def update_coverage_figure(selected_species, rows, selected_rows, plot_species_source):
    if rows == [{}] or rows == [] or rows == None:
        return {"data":[]}
    plot_values = global_vars.plot_values
    traces = []
    trace_ranges = []

     # part of the HACK, replace with "properties.detected_species" when issue is solved
    species_col = "properties_detected_species"
    if plot_species_source == "provided":
        species_col = "properties_provided_species"
    elif plot_species_source == "detected":
        species_col = "properties_detected_species"
    # end HACK

    plot_df = pd.DataFrame(rows)
    if selected_rows is not None and len(selected_rows) > 0:
        plot_df = plot_df.iloc[selected_rows]
    plot_df.loc[pd.isnull(plot_df[species_col]),
                species_col] = "Not classified"

    df_ids = plot_df["_id"]
    if species_col in plot_df.columns and (selected_species in plot_df[species_col].unique() or
        selected_species == "All species"):
        for plot_value in plot_values:
            plot_id = plot_value["id"].replace(".", "_").replace(":", "_")  #HACK
            if selected_species == "All species":
                species_df = plot_df
            else:
                species_df = plot_df[plot_df[species_col] == selected_species]
            species_df[plot_id] = pd.to_numeric(species_df[plot_id], errors="coerce")
            if (plot_id in species_df.columns):
                data_range = plot_value["limits"][1] - plot_value["limits"][0]
                low_limit = min(float(species_df[plot_id].min()), plot_value["limits"][0])
                if low_limit == float(species_df[plot_id].min()):
                    low_limit -= data_range * 0.1
                high_limit = max(float(species_df[plot_id].max()), plot_value["limits"][1])
                if high_limit == float(species_df[plot_id].max()):
                    high_limit += data_range*0.1
                trace_ranges.append([low_limit, high_limit])
                traces.append(
                    go.Box(
                        x=species_df.loc[:, plot_id],
                        text=species_df["name"],
                        marker=dict(
                            size=4
                        ),
                        
                        boxpoints="all",
                        jitter=0.3,
                        pointpos=-1.6,
                        selectedpoints=list(
                            range(len(species_df.index))),
                        name=plot_value["name"],
                        showlegend=False,
                        customdata=species_df["_id"]
                    )
                )
    fig = tools.make_subplots(rows=7, cols=1, print_grid=False)
    fig["layout"].update(
        hovermode="closest",
        title=selected_species,
        height=750,
        margin=go.layout.Margin(
            l=175,
            r=50,
            b=25,
            t=50
        ),
    )
    try:
        fig.append_trace(traces[0], 1, 1)
        fig.append_trace(traces[1], 1, 1)
        fig.append_trace(traces[2], 2, 1)
        fig.append_trace(traces[3], 3, 1)
        fig.append_trace(traces[4], 4, 1)
        fig.append_trace(traces[5], 5, 1)
        fig.append_trace(traces[6], 6, 1)
        fig.append_trace(traces[7], 7, 1)
        
        fig["layout"]["xaxis"].update(range=trace_ranges[0])
        fig["layout"]["xaxis2"].update(range=trace_ranges[2])
        fig["layout"]["xaxis3"].update(range=trace_ranges[3])
        fig["layout"]["xaxis4"].update(range=trace_ranges[4])
        fig["layout"]["xaxis5"].update(range=trace_ranges[5])
        fig["layout"]["xaxis6"].update(range=trace_ranges[6])
        fig["layout"]["xaxis7"].update(range=trace_ranges[7])
    except IndexError:
        pass # It won't draw al ranges/traces because they don't exist.


    fig["layout"]["yaxis"].update(domain=(0.78,1))
    fig["layout"]["yaxis2"].update(domain=(0.655, 0.75))
    fig["layout"]["yaxis3"].update(domain=(0.53, 0.625))
    fig["layout"]["yaxis4"].update(domain=(0.415, 0.5))
    fig["layout"]["yaxis5"].update(domain=(0.28, 0.375))
    fig["layout"]["yaxis6"].update(domain=(0.155, 0.25))
    fig["layout"]["yaxis7"].update(domain=(0.03, 0.125))

    species_size = import_data.get_species_QC_values(selected_species)
    if species_size is None:
        species_size = import_data.get_species_QC_values("default")

    annotations = [
        {
            "x": species_size["min_length"],
            "y": 0,
            "xref": "x",
            "yref": "y",
            "text": "min",
            "arrowhead": 0,
            "ax": 0,
            "ay": 40
        },
        {
            "x": species_size["min_length"],
            "y": 1,
            "xref": "x",
            "yref": "y",
            "text": "min",
            "arrowhead": 0,
            "ax": 0,
            "ay": 40
        },
        {
            "x": species_size["max_length"],
            "y": 0,
            "xref": "x",
            "yref": "y",
            "text": "max",
            "arrowhead": 0,
            "ax": 0,
            "ay": 40
        },
        {
            "x": species_size["max_length"],
            "y": 1,
            "xref": "x",
            "yref": "y",
            "text": "max",
            "arrowhead": 0,
            "ax": 0,
            "ay": 40
        },
        {
            "x": 250000,
            "y": 0,
            "xref": "x2",
            "yref": "y2",
            "text": "max",
            "arrowhead": 0,
            "ax": 0,
            "ay": 40
        },
        { # Cov
            "x": 10,
            "y": 0,
            "xref": "x3",
            "yref": "y3",
            "text": "fail",
            "arrowhead": 0,
            "ax": 0,
            "ay": 40
        },
        { # Cov
            "x": 25,
            "y": 0,
            "xref": "x3",
            "yref": "y3",
            "text": "low",
            "arrowhead": 0,
            "ax": 0,
            "ay": 40
        },
        { # Cov
            "x": 50,
            "y": 0,
            "xref": "x3",
            "yref": "y3",
            "text": "warn",
            "arrowhead": 0,
            "ax": 0,
            "ay": 40
        },
        { # Num reads
            "x": 10000,
            "y": 0,
            "xref": "x5",
            "yref": "y5",
            "text": "min",
            "arrowhead": 0,
            "ax": 0,
            "ay": 40
        },
        {  # Main+uncl
            "x": 0.95,
            "y": 0,
            "xref": "x6",
            "yref": "y6",
            "text": "min",
            "arrowhead": 0,
            "ax": 0,
            "ay": 40
        },
        {  # Uncl
            "x": 0.2,
            "y": 0,
            "xref": "x7",
            "yref": "y7",
            "text": "max",
            "arrowhead": 0,
            "ax": 0,
            "ay": 40
        },
    ]

    fig["layout"].update(annotations=annotations)
    return fig


def create_stamp(value, user):
    return {
        "name": "supplying_lab_check",
        "user-ip": str(request.remote_addr),
        "user": user,
        "date": datetime.datetime.utcnow(),
        "value": value
    }


@app.callback(Output('qc-confirm', 'displayed'),
              [Input('feedback-button', 'n_clicks_timestamp')])
def display_confirm_feedback(button):
    if button is not None:
        return True
    return False


@app.callback(
    Output("placeholder0", "children"),
    [Input("qc-confirm", "submit_n_clicks")],
    [State("qc-user-1", "value")] + [State("sample-radio-{}".format(n), "value")
                                     for n in range(PAGESIZE)]
)
def print_radio(n_clicks_timestamp, user, *args):
    if ("REPORTER_ADMIN" in os.environ and os.environ["REPORTER_ADMIN"] == "True"):
        stamp_a = create_stamp("pass:accepted", user)
        stamp_r = create_stamp("fail:resequence", user)
        stamp_o = create_stamp("fail:other", user)
        stamplist = []
        for val in args:
            if val != "noaction":
                if val.startswith("A_"):
                    stamplist.append((val[2:], stamp_a))
                elif val.startswith("R_"):
                    stamplist.append((val[2:], stamp_r))
                elif val.startswith("O_"):
                    stamplist.append((val[2:], stamp_o))
        if len(stamplist) > 0:
            import_data.email_stamps(stamplist)
            import_data.post_stamps(stamplist)
    return []


server = app.server # Required for gunicorn

if __name__ == '__main__':
    # 0.0.0.0 exposes the app to the network.
    app.run_server(debug=True, host="0.0.0.0", dev_tools_hot_reload=True)
