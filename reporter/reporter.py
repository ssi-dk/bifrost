# -*- coding: utf-8 -*-
import os
import sys
import urllib
from datetime import datetime
from io import StringIO

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import dash_auth
import pandas as pd
import numpy as np
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State

import dash_scroll_up

import flask  #used for the image server

import keys

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
    if species is None or pd.isna(species):
        return None
    words = species.split(" ")
    if len(words) == 1:
        return species
    return "{}. {}".format(words[0][0], " ".join(words[1:]))




app = dash.Dash()

if hasattr(keys, 'USERNAME_PASSWORD'):
    auth = dash_auth.BasicAuth(
        app,
        keys.USERNAME_PASSWORD
    )
app.title = "bifrost"
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
        html.Div(dash_table.DataTable(editable=False), style={"display": "none"}),
        html.H1("bifrost REPORT"),
        html.H2("Loading...", id="run-name"),
        html.Div(id="report-link"),
        dcc.Store(id="data-store", storage_type="memory"),
        html.H4(html.A("Wiki is here", href="https://teams.microsoft.com/l/channel/19%3a7b0b9a088602419e9f84630bacc84c2e%40thread.skype/tab%3a%3a9098abb1-75f5-410a-9011-87db7d42f3c2?label=Wiki&groupId=16852743-838a-400e-921d-6c50cc495b2f&tenantId=d0155445-8a4c-4780-9c13-33c78f22890e")),
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
                        # html.Div(
                        #     [
                        #         html.Button(
                        #             "QCQuickie",
                        #             id="update-qcquickie",
                        #             n_clicks_timestamp=0,
                        #             className="button-primary u-full-width"
                        #         )
                        #     ],
                        #     className="three columns"
                        # ),
                        html.Div(
                            [
                                html.Button(
                                    "Assemblatron",
                                    id="update-assemblatron",
                                    n_clicks_timestamp=0,
                                    className="button-primary u-full-width"
                                )
                            ],
                            className="four columns"
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
                            className="four columns"
                        ),
                        # html.Div(
                        #     [
                        #         html.Button(
                        #             "Table report",
                        #             id="update-table",
                        #             n_clicks_timestamp=0,
                        #             className="button-primary u-full-width"
                        #         )
                        #     ],
                        #     className="three columns"
                        # )
                        html.Div(
                            [
                                html.Button(
                                    "Sample folder",
                                    id="generate-folder",
                                    n_clicks_timestamp=0,
                                    className="button-primary u-full-width"
                                )
                            ],
                            className="four columns"
                        )
                    ],
                    className="row",
                    style={"marginBottom": "15px"}
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
    Output("run-link", "href"),
    [Input("run-list", "value")]
)
def update_run_button(run):
    if run is not None:
        return "/" + run
    else:
        return "/"

@app.callback(
    Output("report-link", "children"),
    [Input("run-name", "children")]
)
def update_run_name(run_name):
    if run_name == "" or run_name == "Not found":
        return None
    else:
        return html.H4(html.A("Link to Run Checker", href="{}/{}".format(keys.run_checker_url, run_name)))

@app.callback(
    Output("report-count", "children"),
    [Input("lasso-sample-ids", "children")]
)
def display_selected_data(ids):
    if ids is not None and len(ids):
        return len(ids.split(","))
    else:
        return 0

@app.callback(
    Output("lasso-div", "children"),
    [Input("summary-plot", "selectedData"),
     Input('datatable-testomatic', 'derived_virtual_data'),
     Input('datatable-testomatic', 'derived_virtual_selected_rows')]
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
        points = [str(point) for point in points if point in union_points]
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
    [Input("run-name", "children"),
    Input("url", "pathname")]
)
def update_group_list(run_name, pathname):
    if run_name == "Loading...":
        return None
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
    if run_name == "Loading...":
        return None
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
        State("data-store", "data")]
        )
def sample_report(page_n, data_content, lasso_selected, data_store):
    if data_content not in ["qcquickie", "assemblatron", "analyzer"] or data_store == "{}":
        return []
    page_n = int(page_n)
    csv_data = StringIO(data_store)
    data = pd.read_csv(csv_data, low_memory=True)
    if lasso_selected != "" and lasso_selected is not None:
        lasso = lasso_selected.split(",")  # lasso first
        data = data[data._id.isin(lasso)]
    if len(data) == 0: return []
    samples = data["_id"]
    #NOTE Could optimize this by not getting all sample's info from mongo before paginating
    data = data.sort_values(["species","name"])
    skips = PAGESIZE * (page_n)
    page = data[skips:skips+PAGESIZE]
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
    [
        # Input("update-qcquickie", "n_clicks_timestamp"),
        Input("update-assemblatron", "n_clicks_timestamp"),
        Input("update-analyzer", "n_clicks_timestamp"),
        # Input("update-table", "n_clicks_timestamp"),
        Input("generate-folder", "n_clicks_timestamp")],
    [State("lasso-sample-ids", "children"),
        State("data-store", "data")]
)
# def update_report(n_qcquickie_ts, n_assemblatron_ts,
def update_report(n_assemblatron_ts,
                n_analyzer_ts, n_generate_ts,
                lasso_selected, data_store):
    n_qcquickie_ts = -1
    n_table_ts = -1
    if lasso_selected != "":
        samples = lasso_selected.split(",")  # lasso first
    elif data_store != None:
        csv_data = StringIO(data_store)
        samples = pd.read_csv(csv_data, low_memory=True)["_id"]
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
    elif n_generate_ts == last_module_ts:
        return generate_sample_folder(samples)
    return []

@app.callback(
    Output(component_id="data-store", component_property="data"),
    [Input(component_id="apply-filter-button", component_property="n_clicks")],
    [State(component_id="species-list", component_property="value"),
        State(component_id="group-list", component_property="value"),
        State(component_id="qc-list", component_property="value"),
        State(component_id="run-name", component_property="children")]
)
def update_selected_samples(n_clicks_ignored, species_list, group_list, qc_list, run_name):
    if run_name == "Loading..." or \
        None in (species_list, group_list, qc_list, run_name):
        return '""'
    else:
        samples = import_data.filter_all(species=species_list,
            group=group_list, qc_list=qc_list, run_name=run_name)
    return samples.to_csv()


@app.callback(
    Output(component_id="testomatic-report",
           component_property="children"), 
    [Input(component_id="data-store", component_property="data")]
)
def update_test_table(data_store):
    empty_table = [
        html.H6('No samples loaded. Click "Apply Filter" to load samples.'),
        html.Div([
            html.P(
                'To filter on a string type eq, space and exact text in double quotes: eq "FBI"'),
            html.P(
                'To filter on a number type eq, < or >, space and num(<number here>): > num(500)')
        ]),
        html.Div(dash_table.DataTable(id="datatable-testomatic",
                                      data=[{}]), style={"display": "none"})
    ]
    if data_store == '""':
        return empty_table
    csv_data = StringIO(data_store)
    tests_df = pd.read_csv(csv_data, low_memory=True)
    if len(tests_df) == 0:
        return empty_table
    qc_action = "stamper:ssi_stamp.assemblatron:action"
    if qc_action not in tests_df:
        tests_df[qc_action] = np.nan
 
    if "R1" not in tests_df:
        tests_df["R1"] = np.nan

    #Temporary fix for Undetermined:
    skipped_mask = tests_df.R1.notnull() & tests_df[qc_action].isnull()
    tests_df.loc[skipped_mask, qc_action] = "skipped"
    no_reads_mask = pd.isnull(tests_df["R1"])
    tests_df.loc[no_reads_mask, qc_action] = "core facility (no reads)"
    mask = pd.isnull(tests_df[qc_action])
    tests_df.loc[mask, qc_action] = "core facility"
    slmask = tests_df[qc_action] == "supplying lab"
    tests_df.loc[slmask, qc_action] = "warning: supplying lab"
    

    # Split test columns
    columns = tests_df.columns
    print(columns)
    split_columns = [
        "stamper:ssi_stamp.assemblatron:1x10xsizediff",
        "stamper:ssi_stamp.whats_my_species:minspecies",
        "stamper:ssi_stamp.whats_my_species:nosubmitted",
        "stamper:ssi_stamp.whats_my_species:detectedspeciesmismatch"
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

    test_cols = [col for col in columns if col.startswith("stamper:ssi_stamp")]
    def concatenate_failed(row):
        res = []
        for col in test_cols:
            test_name = col.split(":")[-1]
            if type(row[col]) == str:
                fields = row[col].split(":")
                if fields[0] in ["fail", "undefined"]:
                    res.append("Test {}: {}, {}".format(test_name, fields[0], fields[1]))
        row["testomatic_failed_tests"] = ". ".join(res)
        return row
    
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
            'width': '250px',
            'padding': '0 15px'
        },
        style_cell_conditional=[
            {
                "if": {"column_id": "testomatic_failed_tests"},
                "textAlign": "left"
            }
        ],
        n_fixed_rows=1,
        row_selectable="multi",
        filtering=True, #Front end filtering
        sorting=True,
        selected_rows=[],
        style_data_conditional=style_data_conditional,
        id="datatable-testomatic"
    )

    rename_dict = {item["id"]:item["name"] for item in COLUMNS}

    renamed = tests_df.rename(rename_dict, axis='columns')#[

    missing_columns = [a for a in list(rename_dict.values()) if not a in list(renamed.columns)]

    # add missing columns
    for column in missing_columns:
        renamed[column] = np.nan
    
    # reorder columns
    renamed = renamed[list(rename_dict.values())]

    csv_string = renamed.to_csv(index=False, encoding="utf-8", sep="\t")
    csv_string = 'data:text/tab-separated-values;charset=utf-8,' + \
        urllib.parse.quote(csv_string)
    return [
        html.H6("Filtered samples ({}):".format(len(tests_df["_id"]))),
        html.Div([
            html.P('To filter on a string type eq, space and exact text in double quotes: eq "FBI"'),
            html.P('To filter on a number type eq, < or >, space and num(<number here>): > num(500)'),
            html.A("Download Table (tsv)", href=csv_string, download='report.tsv')
        ]),
        html.Div(table)
        ]


@app.callback(
    Output(component_id="summary-plot", component_property="selectedData"),
    [Input("data-store", "data"),
     Input(component_id="plot-list", component_property="value"),
     Input('datatable-testomatic', 'derived_virtual_data'),
     Input('datatable-testomatic', 'derived_virtual_selected_rows')]
)
def reset_selection(sample_ids, plot_value, rows, selected_rows):
    return {"points":[]}


@app.callback(
    Output(component_id="summary-plot", component_property="figure"),
    [Input(component_id="plot-list", component_property="value"),
     Input('datatable-testomatic', 'derived_virtual_data'),
     Input('datatable-testomatic', 'derived_virtual_selected_rows')]
)
def update_coverage_figure(plot_value, rows, selected_rows):
    if rows == [{}] or rows == [] or rows == None:
        return {"data":[]}
    plot_query = global_vars.PLOTS[plot_value]["projection"]
    data = []
    plot_df = pd.DataFrame(rows)
    if selected_rows is not None and len(selected_rows) > 0:
        plot_df = plot_df.iloc[selected_rows]
    df_ids = plot_df["_id"]

    species_count = 0
    
    # part of the HACK, replace with "properties.detected_species" when issue is solved
    species_col = "properties_detected_species"
    plot_query = plot_query.replace(".", "_")
    # end HACK

    if species_col in plot_df:

        # reverse the list so it looks right on plot
        species_list = plot_df[species_col].unique()

        for species in reversed(species_list):
            species_df = plot_df[plot_df[species_col] == species]
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
                        selectedpoints=list(
                            range(species_df["_id"].count())),
                        name="{} ({})".format(species_name,
                                              species_df["_id"].count()),
                        showlegend=False,
                        customdata=species_df["_id"]
                    )
            )
    height = max(450, species_count*20 + 200)
    return {
        "data": data,
        "layout": go.Layout(
            hovermode="closest",
            title="{} - selected samples ({})".format(plot_value.replace("_",
                                                               " "), len(plot_df["_id"])),
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
    if run_name == "Loading...":
        return None
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
    if run_name == "Loading...":
        return None
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


@app.callback(
    Output("qc-list", "value"),
    [Input("qc-all", "n_clicks")]
)
def all_QCs(n_clicks):
    return ["OK", "core facility", "supplying lab", "skipped", "Not checked"]

application = app.server # Required for uwsgi

