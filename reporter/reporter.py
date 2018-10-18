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
    if species is None:
        return None
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
        html.H2("Loading...", id="run-name"),
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
    Output("run-selector", "children"),
    [Input("run-name", "children")]
)
def update_run_list(run_name):
    if len(run_name) == 0 or run_name == "Loading...":
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
     Input('datatable-testomatic', 'rows'),
     Input('datatable-testomatic', 'selected_row_indices')]
)
def display_selected_data(selected_data, rows, selected_rows):
    # ignore_this is there so the function is called 
    # when the sample list is updated.
    if rows == [{}] or rows == []:
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
    points = list(dtdf["name"])
    sample_ids = list(dtdf["DB ID"])

    if selected_data is not None and len(selected_data["points"]):
        lasso_points = set([sample["text"]
                    for sample in selected_data["points"]])
        lasso_sample_ids = set([sample["customdata"]
                        for sample in selected_data["points"]])
        union_points = set(points).intersection(lasso_points)
        union_sample_ids = set(sample_ids).intersection(lasso_sample_ids)
        # This way we keep the table order.
        points = [point for point in points if point in union_points]
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
    Output("qc-div", "children"),
    [Input("run-name", "children"),
     Input(component_id="species-list", component_property="value"),
     Input(component_id="group-list", component_property="value")]
)
def update_qc_list(run_name, species, group):
    if run_name == "Loading..." or \
        None in (species, group):
        return dcc.Dropdown(
            id="qc-list",
            multi=True
        )
    num_samples = len(import_data.filter_name(species, group, run_name=run_name))
    if len(run_name) == 0:
        qc_list = import_data.get_qc_list()
    else:
        qc_list = import_data.get_qc_list(run_name)
    qc_options = []
    qc_list_options = []
    sum_items = 0
    for item in qc_list:
        if item["_id"] == None:
            sum_items += item["count"]
            qc_options.append("Not determined")
            qc_list_options.append({
                "label": "Not determined ({})".format(item["count"]),
                "value": "Not determined"
            })
        else:
            sum_items += item["count"]
            qc_options.append(item["_id"])
            qc_list_options.append({
                "label": "{} ({})".format(item["_id"], item["count"]),
                "value": item["_id"]
            })
    if sum_items < num_samples:
            qc_options.append("Not tested")
            qc_list_options.append({
                "label": "Not tested ({})".format(num_samples - sum_items),
                "value": "Not tested"
            })
    return dcc.Dropdown(
        id="qc-list",
        options=qc_list_options,
        multi=True,
        value=qc_options
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
        State("selected-samples-ids", "children")]
        )
def sample_report(page_n, data_content, lasso_selected, prefilter_samples):
    if data_content not in ["qcquickie", "assemblatron", "analyzer"]:
        return []
    page_n = int(page_n)
    if lasso_selected != "" and lasso_selected is not None:
        samples = lasso_selected.split(",")  # lasso first
    else:
        samples = prefilter_samples.split(",")
    if samples == [""]: return []
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
    [
        # Input("update-qcquickie", "n_clicks_timestamp"),
        Input("update-assemblatron", "n_clicks_timestamp"),
        Input("update-analyzer", "n_clicks_timestamp"),
        # Input("update-table", "n_clicks_timestamp"),
        Input("generate-folder", "n_clicks_timestamp")],
    [State("lasso-sample-ids", "children"),
        State("selected-samples-ids", "children")]
)
# def update_report(n_qcquickie_ts, n_assemblatron_ts,
def update_report(n_assemblatron_ts,
                n_analyzer_ts, n_generate_ts,
                lasso_selected, prefilter_samples):
    n_qcquickie_ts = -1
    n_table_ts = -1
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
        Input(component_id="qc-list", component_property="value"),
        Input(component_id="run-name", component_property="children")]
)
def update_selected_samples(species_list, group_list, qc_list, run_name):
    if run_name == "Loading..." or \
        None in (species_list, group_list, qc_list, run_name):
        return None
    samples = import_data.filter_name(species=species_list,
        group=group_list, qc_list=qc_list, run_name=run_name)
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
    Output(component_id="testomatic-report",
           component_property="children"), 
    [Input(component_id="selected-samples-list", component_property="value")],
     [State(component_id="species-list", component_property="value"),
        State(component_id="group-list", component_property="value"),
        State(component_id="run-name", component_property="children"),
        State(component_id="qc-list", component_property="value")]
)
def update_test_table(selected_samples, species_list, group_list, run_name, qc_list):
    if run_name == "Loading..." or run_name == "" or \
            None in (species_list, group_list, qc_list, run_name):
        #return None
        return [html.H6("Filtered samples (0):"), dt.DataTable(id="datatable-testomatic", rows=[{}])]
    columns = ["name", 'testomatic.assemblatron:action',  "sample_sheet.Comments",
            "sample_sheet.group", "properties.provided_species", "properties.detected_species",
            'testomatic.assemblatron:1xgenomesize',
            'testomatic.assemblatron:10xgenomesize',
            'testomatic.assemblatron:1x10xsizediff', 'testomatic.assemblatron:avgcoverage',
            'assemblatron.bin_contigs_at_1x',
            'assemblatron.snp_filter_10x_10%',
            'testomatic.assemblatron:numreads',
            # 'testomatic.qcquickie:action',
            # 'testomatic.qcquickie:1xgenomesize',
            # 'testomatic.qcquickie:10xgenomesize',
            # 'testomatic.qcquickie:1x10xsizediff',
            # 'testomatic.qcquickie:avgcoverage',
            # 'qcquickie.bin_contigs_at_1x',
            # 'testomatic.qcquickie:numreads',
            'testomatic.whats_my_species:maxunclassified',
            'testomatic.whats_my_species:minspecies',
            'testomatic.whats_my_species:nosubmitted',
            'testomatic.whats_my_species:submitted==detected', "_id"]
    column_names = ['name', 'QC action', "Comments", 'Supplying lab', 'Provided Species',
                    'Detected Species', 'Genome size (1x)', 'Genome size (10x)',
                    'G. size difference (1x - 10x)',
                    'Avg. coverage', '# contigs',
                    'Ambiguous sites',
                    '# reads',
                    # 'qcquickie QC',
                    # 'qcquickie 1xgenomesize', 'qcquickie 10xgenomesize',
                    # 'qcquickie 1x10xsizediff', 'qcquickie avgcoverage',
                    # 'qcquickie # contigs',
                    # 'qcquickie numreads',
                    'Unclassified reads',
                    'Main species read %', 'Detected species in DB',
                    'Submitted sp. is same as detected', 'DB ID']
    tests_df = import_data.filter_all(
        species_list, group_list, qc_list, run_name)
    
    tests_df = tests_df.reindex(columns=columns)
    tests_df.columns = column_names
    mask = pd.isnull(tests_df["QC action"])
    tests_df.loc[mask, "QC action"] = "core facility"
    slmask = tests_df["QC action"] == "supplying lab"
    tests_df.loc[slmask, "QC action"] = "warning: supplying lab"

    table = dt.DataTable(
        rows=tests_df.to_dict("records"),

        # columns=global_vars.columns, # sets the order
        column_widths=[175] * len(columns),
        row_selectable=True,
        editable=False,
        filterable=True,
        sortable=True,
        selected_row_indices=[],
        id="datatable-testomatic"
    )

    # th = html.Thead(
    #     html.Tr(list(map(lambda x: html.Th(x), column_names)), className="trow header"))
    # tbody = []
    # for sample in tests_df.iterrows():
    #     row = []
    #     for value, column in zip(sample[1][columns], columns):
    #         if str(value).startswith("fail") or str(value).startswith("undefined") \
    #         or value == "supplying lab" or value =="core facility":
    #             td = html.Td(str(value), className="cell red")
    #         elif str(value).startswith("KeyError") or (pd.isnull(value) and column.endswith("action")):
    #             td = html.Td(str(value), className="cell yellow")
    #         else:
    #             td = html.Td(str(value), className="cell")
    #         row.append(td)
    #     tbody.append(html.Tr(row, className="trow"))
    # tb = html.Tbody(tbody)

    # tests_df = tests_df[columns]
    
    #tests_df.columns = column_names
    csv_string = tests_df.to_csv(index=False, encoding="utf-8", sep="\t")
    csv_string = 'data:text/tab-separated-values;charset=utf-8,' + \
        urllib.parse.quote(csv_string)
    return [
        html.H6("Filtered samples ({}):".format(len(tests_df["DB ID"]))),
        html.A("Download Table (tsv)", href=csv_string, download='report.tsv'),
        table
        # html.Table([th, tb], className="fixed-header")
        ]


@app.callback(
    Output(component_id="summary-plot", component_property="selectedData"),
    [Input("selected-samples-ids", "children"),
     Input(component_id="plot-list", component_property="value"),
     Input('datatable-testomatic', 'rows'),
     Input('datatable-testomatic', 'selected_row_indices')]
)
def reset_selection(sample_ids, plot_value, rows, selected_rows):
    return {"points":[]}


@app.callback(
    Output(component_id="summary-plot", component_property="figure"),
    [Input("selected-samples-ids", "children"),
     Input(component_id="plot-list", component_property="value"),
     Input('datatable-testomatic', 'rows'),
     Input('datatable-testomatic', 'selected_row_indices')]
)
def update_coverage_figure(sample_ids, plot_value, rows, selected_rows):
    samples = sample_ids.split(",")
    if rows == [{}] or rows == [] or (len(samples) == 1 and samples[0] == ""):
        return {"data":[]}
    plot_query = global_vars.PLOTS[plot_value]["projection"]
    plot_func = global_vars.PLOTS[plot_value].get("func")

    data = []
    plot_df = import_data.filter_all(
        sample_ids=sample_ids.split(","), func=plot_func)

    dtdf = pd.DataFrame(rows)
    if selected_rows is not None and len(selected_rows) > 0:
        dtdf = dtdf.iloc[selected_rows]
        df_ids = dtdf["_id"]
        plot_df = plot_df[plot_df._id.isin(df_ids)]

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
                                                               " "), len(dtdf["DB ID"])),
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
    [Input("qc-all", "n_clicks")],
    [State("run-name", "children")]
)
def all_QCs(n_clicks, run_name):
    if run_name == "Loading...":
        return None
    if len(run_name) == 0:
        qc_list = import_data.get_qc_list()
    else:
        qc_list = import_data.get_qc_list(run_name)
    qc_options = []
    for item in qc_list:
        if item["_id"] == None:
            qc_options.append("Not determined")
        else:
            qc_options.append(item["_id"])
    qc_options.append("Not tested")
    return qc_options

application = app.server # Required for uwsgi

