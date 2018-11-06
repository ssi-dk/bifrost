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
        html.Div(dash_table.DataTable(editable=False), style={"display": "none"}),
        html.H1("SerumQC REPORT"),
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
     Input('datatable-testomatic', 'selected_rows')]
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
    sample_ids = list(dtdf["DB_ID"])

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


# @app.callback(
#     Output("qc-div", "children"),
#     [Input("run-name", "children"),
#      Input(component_id="species-list", component_property="value"),
#      Input(component_id="group-list", component_property="value")]
# )
# def update_qc_list(run_name, species, group):
#     if run_name == "Loading..." or \
#         None in (species, group):
#         return dcc.Dropdown(
#             id="qc-list",
#             multi=True
#         )
#     num_samples = len(import_data.filter_name(species, group, run_name=run_name))
#     if len(run_name) == 0:
#         qc_list = import_data.get_qc_list()
#     else:
#         qc_list = import_data.get_qc_list(run_name)
#     qc_options = []
#     qc_list_options = []
#     sum_items = 0
#     for item in qc_list:
#         if item["_id"] == None:
#             sum_items += item["count"]
#             qc_options.append("Not determined")
#             qc_list_options.append({
#                 "label": "Not determined ({})".format(item["count"]),
#                 "value": "Not determined"
#             })
#         else:
#             sum_items += item["count"]
#             qc_options.append(item["_id"])
#             qc_list_options.append({
#                 "label": "{} ({})".format(item["_id"], item["count"]),
#                 "value": item["_id"]
#             })
#     if sum_items < num_samples:
#             qc_options.append("Not tested")
#             qc_list_options.append({
#                 "label": "Not tested ({})".format(num_samples - sum_items),
#                 "value": "Not tested"
#             })
#     return dcc.Dropdown(
#         id="qc-list",
#         options=qc_list_options,
#         multi=True,
#         value=qc_options
#     )


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
        return "{}"
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
    columns = ['testomatic.assemblatron:action', "name",  "sample_sheet.Comments",
            "sample_sheet.group", "properties.provided_species", "properties.detected_species",
            'testomatic.assemblatron:1xgenomesize',
            'testomatic.assemblatron:10xgenomesize',
            'testomatic.assemblatron:1x10xsizediff', 'testomatic.assemblatron:avgcoverage',
            'assemblatron.bin_contigs_at_1x',
            'assemblatron.snp_filter_10x_10%',
            'testomatic.assemblatron:numreads',
            "analyzer.mlst_report",
            'testomatic.whats_my_species:maxunclassified',
            'testomatic.whats_my_species:minspecies',
            'testomatic.whats_my_species:nosubmitted',
            'testomatic.whats_my_species:submitted==detected', "_id"]
    column_names = ['QC_action', 'name', "Comments", 'Supplying_lab', 'Provided_Species',
                    'Detected_Species', 'Genome_size_1x', 'Genome_size_10x',
                    'G_size_difference_1x_10',
                    'Avg._coverage', 'num_contigs',
                    'Ambiguous_sites',
                    'num_reads',
                    'mlst',
                    'Unclassified_reads',
                    'Main_species_read_percent', 'Detected_species_in_DB',
                    'Submitted_sp_is_same_as_detected', 'DB_ID']
    if data_store == "{}":
        return [
            html.H6("Filtered samples (0):"),
            dash_table.DataTable(
                data=[],
                style_table={
                    'overflowX': 'scroll',
                    'overflowY': 'scroll',
                    'maxHeight': '480'
                },
                columns=[{"name": i, "id": i} for i in columns],
                # n_fixed_columns=1,
                style_cell={'width': '150px'},

                # n_fixed_rows=1, # NOT WORKING NOTE
                row_selectable="multi",
                filtering=True,  # Front end filtering
                sorting=True,
                selected_rows=[],
                id="datatable-testomatic"
            )
        ]
    csv_data = StringIO(data_store)
    tests_df = pd.read_csv(csv_data, low_memory=True)
    tests_df = tests_df.reindex(columns=columns)
    tests_df.columns = column_names
    mask = pd.isnull(tests_df["QC_action"])
    tests_df.loc[mask, "QC_action"] = "core facility"
    slmask = tests_df["QC_action"] == "supplying lab"
    tests_df.loc[slmask, "QC_action"] = "warning: supplying lab"
    print(tests_df.to_dict("rows"))
    table = dash_table.DataTable(

        data=tests_df.to_dict("rows"),
        style_table={
            'overflowX': 'scroll',
            'overflowY': 'scroll',
            'maxHeight': '480'
        },
        columns=[{"name": i, "id": i} for i in tests_df.columns],
        # n_fixed_columns=1,
        style_cell={'width': '150px'},
        
        # n_fixed_rows=1, # NOT WORKING NOTE
        row_selectable="multi",
        filtering=True, #Front end filtering
        sorting=True,
        selected_rows=[],
        id="datatable-testomatic"
    )

    csv_string = tests_df.to_csv(index=False, encoding="utf-8", sep="\t")
    csv_string = 'data:text/tab-separated-values;charset=utf-8,' + \
        urllib.parse.quote(csv_string)
    return [
        html.H6("Filtered samples ({}):".format(len(tests_df["DB_ID"]))),
        html.A("Download Table (tsv)", href=csv_string, download='report.tsv'),
        table
        # html.Table([th, tb], className="fixed-header")
        ]


@app.callback(
    Output(component_id="summary-plot", component_property="selectedData"),
    [Input("data-store", "data"),
     Input(component_id="plot-list", component_property="value"),
     Input('datatable-testomatic', 'derived_virtual_data'),
     Input('datatable-testomatic', 'selected_rows')]
)
def reset_selection(sample_ids, plot_value, rows, selected_rows):
    return {"points":[]}


@app.callback(
    Output(component_id="summary-plot", component_property="figure"),
    [Input(component_id="plot-list", component_property="value"),
     Input('datatable-testomatic', 'derived_virtual_data'),
     Input('datatable-testomatic', 'selected_rows')],
    [State('data-store', 'data')]
)
def update_coverage_figure(plot_value, rows, selected_rows, data_store):
    if rows == [{}] or rows == [] or rows == None:
        return {"data":[]}
    plot_query = global_vars.PLOTS[plot_value]["projection"]

    data = []
    csv_data = StringIO(data_store)
    plot_df = pd.read_csv(csv_data, low_memory=True)
    if selected_rows is not None and len(selected_rows) > 0:
        plot_df = plot_df.iloc[selected_rows]
    df_ids = plot_df["_id"]

    species_count = 0
    if 'species' in plot_df:

        # reverse the list so it looks right on plot
        species_list = plot_df.species.unique()

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
                        selectedpoints=list(range(species_df["_id"].count())),
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
    return ["OK", "core facility", "supplying lab", "Not tested"]

application = app.server # Required for uwsgi

