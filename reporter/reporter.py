# -*- coding: utf-8 -*-
import os
import sys
import datetime

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dt
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State

import flask  #used for the image server
import glob

import mongo_interface

#DEV images
image_directory = '/Users/mbas/Documents/SerumQC-private/reporter/resources/img/'
list_of_images = [os.path.basename(x) for x in glob.glob(
    '{}*.svg'.format(image_directory))]
static_image_route = '/static/'


def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

# Globals
COLOR_DICT = mongo_interface.get_species_colors()

PLOT_VALUES = [
    "bin_length_at_1x",
    "bin_length_at_10x",
    "bin_length_at_25x",
    "bin_length_1x_25x_diff",
    "bin_contigs_at_1x",
    "bin_contigs_at_10x",
    "bin_contigs_at_25x",
    "N50",
    "N75"
]
DEFAULT_PLOT = 0

PAGESIZE = 50

def get_species_color(species):
    color = COLOR_DICT.get(species, "#b3ccc1")
    if color == '':
        color = "#b3ccc1"  # Default
    return color

def filter_dataframe(dataframe, species_list, group_list, page_n=None):
    filtered = dataframe[
        dataframe.short_class_species_1.isin(species_list) &
        dataframe.supplying_lab.isin(group_list)
    ]
    if page_n is None:
        return filtered
    else:
        return filtered.iloc[page_n*PAGESIZE:(page_n+1)*PAGESIZE]


def format_selected_samples(filtered_df):
    "Returns a formatted string of selected samples"
    return "\n".join([row["name"] for index, row in filtered_df.iterrows()])


# Data loading before application starts.

def get_qcquickie(sample_db):
    sample = {}
    sample['_id'] = str(sample_db['_id'])
    sample["input_read_status"] = sample_db["sample"]["input_read_status"]
    sample['name'] = sample_db['sample'].get('name')
    sample['user'] = sample_db['sample'].get('user')
    sample['R1_location'] = sample_db['sample'].get('R1')
    sample['R2_location'] = sample_db['sample'].get('R2')
    sample["run_name"] = sample_db["sample"]["run_folder"].split("/")[-1]
    if "setup_time" in sample_db["sample"]:
        sample["setup_time"] = datetime.datetime.strptime(sample_db["sample"]["setup_time"], "%Y-%m-%d %H:%M:%S.%f")
    else:
        sample["setup_time"] = None

    # Not nice, but we have many different combinations
    if "sample_sheet" in sample_db["sample"]:
        if sample['name'] == None:
            sample['name'] = sample_db['sample']["sample_sheet"]["sample_name"]
        sample["provided_species"] = sample_db['sample']["sample_sheet"].get("provided_species", "")
        sample['supplying_lab'] = sample_db['sample']["sample_sheet"]['group']
        sample['comments'] = sample_db["sample"]["sample_sheet"]["Comments"]
        sample["emails"] = sample_db["sample"]["sample_sheet"]["emails"]
    

    if "qcquickie" in sample_db:
        for key, value in sample_db['qcquickie']['summary'].items():
            sample[key] = value
        sample["N50"] = sample_db["qcquickie"]["quast/report_tsv"]["N50"]
        sample["N75"] = sample_db["qcquickie"]["quast/report_tsv"]["N75"]
        sample["bin_length_1x_25x_diff"] = sample["bin_length_at_1x"] - \
            sample["bin_length_at_25x"]
    
    if not "run_name" in sample:
        sys.stderr.write("Sample {} has no run name.\n".format(sample["name"]))
        sample["run_name"] = ""

    if "name_classified_species_1" in sample:
        genus, species = sample["name_classified_species_1"].split()
        sample["short_class_species_1"] = '{}. {}'.format(genus[0], species)
    else:
        sample["name_classified_species_1"] = "Not classified"
        sample["short_class_species_1"] = "Not classified"

    if not "supplying_lab" in sample:
        sample["supplying_lab"] = "Not specified"

    return sample

# Components


def html_table(data, **kwargs):
    return html.Table([
        html.Tr([html.Td(data_cell) for data_cell in data_row])
        for data_row in data
    ], **kwargs)


def html_td_percentage(value, color):
    string = str(round(float(value) * 100, 2)) + "%"
    return html.Td(
        html.Div([
            html.Span(string, className="val"),
            html.Span(
                className="bar",
                style={"backgroundColor": color, "width": string}
            )
        ], className="wrapper"
    ), className="data-colored")


# QCQuickie

def generate_sample_report(dataframe, sample):
    return (
        html.Div(
            [
                html.A(id="sample-" + sample["name"]),
                html.H5(
                    sample["name"],
                    className="box-title"
                ),
                html_sample_tables(sample, className="row"),

                graph_sample_depth_plot(
                    sample, dataframe[dataframe.name_classified_species_1 == sample["name_classified_species_1"]])
            ],
            className="border-box"
        )
    )


def html_species_report(dataframe, species, **kwargs):
    report = []
    for index, sample in dataframe.loc[dataframe["name_classified_species_1"] == species].iterrows():
        report.append(generate_sample_report(dataframe, sample))
    return html.Div(report, **kwargs)


def html_organisms_table(sample_data, **kwargs):
    percentages = [
        sample_data["percent_classified_species_1"],
        sample_data["percent_classified_species_2"],
        sample_data["percent_unclassified"]
    ]

    color_1 = get_species_color(sample_data["name_classified_species_1"])  # Default
#    color_2 = COLOR_DICT.get(
#        sample_data["name_classified_species_2"], "#f3bbd3")  # Default
    color_2 = "#f3bbd3"  # Default

#   color_u = COLOR_DICT.get("", "#fee1cd")  # Default
    color_u = "#fee1cd"  # Default

    return html.Div([
        html.H6("Detected Organisms", className="table-header"),
        html.Table([
            html.Tr([
                html.Td(html.I(sample_data["name_classified_species_1"])),
                html_td_percentage(percentages[0], color_1)
            ]),
            html.Tr([
                html.Td(html.I(sample_data["name_classified_species_2"])),
                html_td_percentage(percentages[1], color_2)
            ]),
            html.Tr([
                html.Td("Unclassified"),
                html_td_percentage(percentages[2], color_u)
            ])
        ])
    ], **kwargs)


def html_sample_tables(sample_data, **kwargs):
    """Generate the tables for each sample containing submitter information,
       detected organisms etc. """
    genus = str(sample_data["name_classified_species_1"]).split()[0].lower()
    if "{}.svg".format(genus) in list_of_images:
        img = html.Img(src="/static/" + str(sample_data["name_classified_species_1"]).split()
                       [0].lower() + ".svg", className="svg_bact")
    else:
        img = []
    if type(sample_data["emails"]) is str:
        n_emails = len(sample_data["emails"].split(";"))
        if (n_emails > 1):
            emails = ", ".join(sample_data["emails"].split(";")[:2])
            if (n_emails > 2):
                emails += ", ..."
        else:
            emails = sample_data["emails"]
    else:
        emails = ''

    return html.Div([
        html.Div(img, className="bact_div"),
        html.H6("Run folder: " + sample_data["run_name"]),
        html.H5("Sample Data", className="table-header"),
        html.Div([
            html.Div([
                html_table([
                    ["Sample name", sample_data["name"]],
                    ["Supplying lab", sample_data["supplying_lab"]],
                    ["Submitter emails", emails],
                    ["Submitter user", sample_data["user"]],
                    ["Setup time", sample_data["setup_time"]],
                    ["Provided species", html.I(
                        sample_data["provided_species"])]
                ])
            ], className="six columns"),
            html.Div([
                html.H6("User Comments", className="table-header"),
                sample_data["comments"]
            ], className="six columns"),
        ], className="row"),
        html.H5("QCQuickie Results", className="table-header"),
        html.Div([
            html.Div([
                html_table([
                    [
                        "Number of contigs",
                        "{:,}".format(sample_data["bin_contigs_at_1x"])
                    ],
                    [
                        "N50",
                        "{:,}".format(sample_data["N50"])
                    ],
                    [
                        "N75",
                        "{:,}".format(sample_data["N75"])
                    ],
                    [
                        "bin length at 1x depth",
                        "{:,}".format(sample_data["bin_length_at_1x"])
                    ],
                    [
                        "bin length at 10x depth",
                        "{:,}".format(sample_data["bin_length_at_10x"])
                    ],
                    [
                        "bin length at 25x depth",
                        "{:,}".format(sample_data["bin_length_at_25x"])
                    ]
                ])
            ], className="six columns"),
            html_organisms_table(sample_data, className="six columns")
        ], className="row")
    ], **kwargs)


def graph_sample_depth_plot(sample, background_dataframe):
    # With real data, we should be getting sample data (where to put 1, 10
    # and 25x annotation) and the info for the rest of that species box.
    return dcc.Graph(
        id="coverage-1-" + sample['name'],
        figure={
            "data": [
                go.Box(
                    x=background_dataframe["bin_length_at_1x"],
                    text=background_dataframe["name"],
                    boxpoints="all",
                    jitter=0.3,
                    pointpos=-1.8,
                    marker=dict(color=get_species_color(
                        sample['name_classified_species_1']))
                )
                #{"x": [1, 2, 3], "y": [2, 4, 5], "type": "bar", "name": u"Montréal"},
            ],
            "layout": go.Layout(
                title="{}: Binned Depth 1x size".format(sample['name']),
                margin=go.Margin(
                    l=75,
                    r=50,
                    b=25,
                    t=50
                ),
                annotations=go.Annotations([
                    go.Annotation(
                        x=sample["bin_length_at_1x"],
                        y=0,
                        text="1x",
                        showarrow=True,
                        ax=40,
                        ay=0
                    ),
                    go.Annotation(
                        x=sample["bin_length_at_10x"],
                        y=0.0,
                        text="10x",
                        showarrow=True,
                        ax=0,
                        ay=40
                    ),
                    go.Annotation(
                        x=sample["bin_length_at_25x"],
                        y=0,
                        text="25x",
                        showarrow=True,
                        ax=0,
                        ay=-40
                    ),
                ])
            )
        },
        style={"height": "200px"}
    )


def children_sample_list_report(filtered_df):
    report = []
    for species in filtered_df.name_classified_species_1.unique():
        report.append(html.Div([
            html.A(id="species-cat-" + str(species).replace(" ", "-")),
            html.H4(html.I(str(species))),
            html_species_report(filtered_df, species)
        ]))
    return report


def html_table_run_information(run_data):
    return html_table([
        ["Run Name", run_data["run_name"]],
        ["Run Date", "06 MAY 2018"],
        ["Placeholder", "Lorem ipsum"]
    ])


def html_div_summary(dataframe, species_list, sample_list, group_list):
    species_list = dataframe.short_class_species_1.unique()
    species_list_options = [{"label": species, "value": species}
                            for species in species_list]

    sample_list = dataframe.name.unique()
    sample_list_options = [{"label": sample, "value": sample}
                           for sample in sample_list]

    group_list = dataframe.supplying_lab.unique()
    group_list_options = [{"label": group, "value": group}
                          for group in group_list]

    plot_values_options = [{"label": plot, "value": plot}
                           for plot in PLOT_VALUES]
    return html.Div(
        [
            html.H5("Summary", className="box-title"),
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            html.Label(
                                                [
                                                    "Supplying lab ",
                                                    html.Small(
                                                        [
                                                            "(",
                                                            html.A(
                                                                "all",
                                                                href="#",
                                                                n_clicks=0,
                                                                id="group-all"
                                                            ),
                                                            ")"
                                                        ]
                                                    )
                                                ],
                                                htmlFor="group-list"
                                            ),
                                            dcc.Dropdown(
                                                id="group-list",
                                                options=group_list_options,
                                                multi=True,
                                                value=group_list
                                            )
                                        ],
                                        className="twelve columns"
                                        )
                                ],
                                className="row"
                            ),
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            html.Label(
                                                [
                                                    "Detected organism ",
                                                    html.Small(
                                                        [
                                                            "(",
                                                            html.A(
                                                                "all",
                                                                href="#",
                                                                n_clicks=0,
                                                                id="species-all"
                                                            ),
                                                            ")"
                                                        ]
                                                    )
                                                ],
                                                htmlFor="group-list"
                                            ),
                                            dcc.Dropdown(
                                                id="species-list",
                                                options=species_list_options,
                                                multi=True,
                                                value=species_list
                                            )
                                        ],
                                        className="twelve columns"
                                    )
                                ],
                                className="row"
                            ),
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            html.Label("Plot value",
                                                       htmlFor="plot-list"),
                                            dcc.Dropdown(
                                                id="plot-list",
                                                options=plot_values_options,
                                                value=PLOT_VALUES[DEFAULT_PLOT]
                                            )
                                        ],
                                        className="twelve columns"
                                    )
                                ],
                                className="row"
                            )
                        ],
                        className="eight columns"
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            html.Label(
                                                [
                                                    "Selected Samples (",
                                                    html.Span(id="sample-count"),
                                                    "):"
                                                ],
                                                       htmlFor="plot-list"),
                                            dcc.Textarea(
                                                id="selected-samples",
                                                className="u-full-width",
                                                style={"resize": "none",
                                                       "height": "400px"},
                                                readOnly=True
                                            )
                                        ],
                                        className="twelve columns"
                                    )
                                ],
                                className="row"
                            )
                        ],
                        className="four columns"
                    )
                ],
                className="row"
            ),
            # Filtered list
            html.Div(
                [
                    html.Div(
                        [
                        ],
                        className="four columns"
                    ),
                    html.Div(
                        [
                            html.Label("Go to sample", htmlFor="sample-list"),
                            html.Div([
                                html.Div(
                                    dcc.Dropdown(
                                        id="sample-list",
                                        options=sample_list_options,
                                        placeholder="Sample name"
                                    )
                                , className="six columns"),
                                html.Div(
                                    html.A("Go", id="go-to-sample",
                                        href="#", className="button")
                                , className="six columns")
                            ], className="row")
                        ],
                        className="six columns"
                    ),
                ], className="row"
            ),

            dcc.Graph(id="summary-plot")
        ], className="border-box"
    )


def main(argv):
    run_folder = "../testdata/" + "tiny_180608"
    dataframe = pd.DataFrame(list(map(get_qcquickie, mongo_interface.test_get_all_samples())))
    

    # Change this to provided species ?
    dataframe = dataframe.sort_values("name_classified_species_1")
    app = dash.Dash()
    app.config['suppress_callback_exceptions'] = True

    species_list = dataframe.short_class_species_1.unique()

    sample_list = dataframe.name.unique()

    group_list = dataframe.supplying_lab.unique()

    # Temp css to make it look nice
    # Dash CSS
    app.css.append_css(
        {"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})
    # Custom CSS, host elsewere for production
    app.css.append_css(
        {"external_url": "https://codepen.io/martinbaste/pen/JZKeRv.css"})
    # Lato font
    app.css.append_css(
        {"external_url": "https://fonts.googleapis.com/css?family=Lato"})

    app.layout = html.Div(className="container", children=[
        html.Div(dt.DataTable(rows=[{}], editable=False), style={'display': 'none'}),
        html.A(id="top"),
        html.H1("QC REPORT"),
        html.H2(run_folder),
        dcc.Location(id="url", refresh=False),
        html_table_run_information({'run_name': run_folder}),
        html_div_summary(dataframe, species_list, sample_list, group_list),
        html.Div(
            [
                html.H5("Update report", className="box-title"),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Button(
                                    "Individual sample report",
                                    id="update-samples",
                                    n_clicks_timestamp=0,
                                    className="button-primary u-full-width"
                                )
                            ],
                            className="six columns"
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
                            className="six columns"
                        )
                    ],
                    className="row"
                )
            ],
            className="border-box"
        ),
        html.Div(id="current-report")
    ])


# Callbacks

    # We could make this one much faster by hiding the unused species with CSS
    # by adding a new hidden class.

    #While dev

    @app.server.route('{}<image_path>.svg'.format(static_image_route))
    def serve_image(image_path):
        image_name = '{}.svg'.format(image_path)
        if image_name not in list_of_images:
            raise Exception(
                '"{}" is excluded from the allowed static files'.format(image_path))
        return flask.send_from_directory(image_directory, image_name)

    @app.callback(
        Output(component_id="page-n",
               component_property="children"),
        [Input(component_id="prevpage", component_property="n_clicks_timestamp"),
         Input(component_id="nextpage", component_property="n_clicks_timestamp")],
        [State(component_id="page-n", component_property="children"),
         State(component_id="species-list", component_property="value"),
         State(component_id="group-list", component_property="value")]
    )
    def next_page(prev_ts, next_ts, page_n, species_list, group_list):
        page_n = int(page_n)
        max_page = filter_dataframe(dataframe, species_list, group_list)["_id"].count() // PAGESIZE
        if prev_ts > next_ts:
            return max(page_n - 1, 0)
        elif next_ts > prev_ts:
            return min(page_n + 1, max_page)
        else:
            return 0

    @app.callback(
        Output(component_id="sample-report", component_property="children"),
        [Input(component_id="page-n", component_property="children")],
         [State(component_id="species-list", component_property="value"),
          State(component_id="group-list", component_property="value")]
         )
    def sample_report(page_n, species_list, group_list):
        page_n = int(page_n)
        filtered = filter_dataframe(dataframe, species_list, group_list, page_n)
        max_page = len(filter_dataframe(dataframe, species_list, group_list)) // PAGESIZE
        return [
            html.H4("Page {} of {}".format(page_n + 1, max_page + 1)),
            html.Div(children_sample_list_report(filtered))
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
        [Input(component_id="page-n", component_property="children")],
        [State(component_id="species-list", component_property="value"),
         State(component_id="group-list", component_property="value")]
    )
    def update_nextpage(page_n, species_list, group_list):
        page_n = int(page_n)
        max_page = len(filter_dataframe(
            dataframe, species_list, group_list)) // PAGESIZE
        if page_n == max_page:
            return True
        else:
            return False

    @app.callback(
        Output(component_id="current-report", component_property="children"),
        [Input(component_id="update-samples", component_property="n_clicks_timestamp"),
         Input(component_id="update-table", component_property="n_clicks_timestamp")]

    )
    def update_report(n_samples_ts, n_table_ts):
        if n_samples_ts > n_table_ts:  # samples was clicked
            return [
                html.H3("Individual Sample Reports"),
                html.Span("0", style={'display': 'none'}, id="page-n"),
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
            ]
        elif n_table_ts > n_samples_ts:  # table was clicked
            columns = ['name', 'supplying_lab', 'run_name', '_id', 'input_read_status',
                       'emails', 'user', 'R1_location', 'R2_location', 'provided_species',
                       'name_classified_species_1', 'percent_classified_species_1',
                       'name_classified_species_2', 'percent_classified_species_2',
                       'percent_unclassified', 'bin_length_at_1x', 'bin_length_at_10x',
                       'bin_length_at_25x', 'bin_length_1x_25x_diff', 'bin_coverage_at_1x',
                       'bin_coverage_at_10x', 'bin_coverage_at_25x', 'bin_contigs_at_1x',
                       'bin_contigs_at_10x', 'N50', 'N75', 'snp_filter_deletions',
                       'snp_filter_indels', 'snp_filter_10x_10%', 'comments']
            return [
                html.H3("Table Report"),
                dt.DataTable(
                    rows=filter_dataframe(
                        dataframe, species_list, group_list).to_dict("records"),

                    # optional - sets the order of columns
                    columns=columns,
                    column_widths=[150]*len(columns),
                    editable=False,
                    filterable=True,
                    sortable=True,
                    selected_row_indices=[],
                    id='datatable-samples'
                )
            ]
        return []

    @app.callback(
        Output(component_id="selected-samples", component_property="value"),
        [Input(component_id="species-list", component_property="value"),
         Input(component_id="group-list", component_property="value")]
    )
    def update_selected_samples(species_list, group_list):
        filtered_df = filter_dataframe(dataframe, species_list, group_list)
        return [format_selected_samples(filtered_df)]

    @app.callback(
        Output(component_id="sample-count",
                component_property="children"),
        [Input(component_id="species-list", component_property="value"),
            Input(component_id="group-list", component_property="value")]
    )
    def update_sample_count(species_list, group_list):
        filtered_df = filter_dataframe(dataframe, species_list, group_list)
        return filtered_df['name'].count()

    @app.callback(
        Output(component_id="go-to-sample", component_property="href"),
        [Input(component_id="sample-list", component_property="value")]
    )
    def update_go_to_sample(sample_name):
        if sample_name in sample_list:
            return "#sample-" + sample_name

    @app.callback(
        Output(component_id="summary-plot", component_property="figure"),
        [Input(component_id="species-list", component_property="value"),
         Input(component_id="group-list", component_property="value"),
         Input(component_id="plot-list", component_property="value")]
    )
    def update_coverage_figure(species_list, group_list, plot_value):
        data = []
        filtered_df = filter_dataframe(dataframe, species_list, group_list)
        for species in species_list:
            species_df = filtered_df[filtered_df.short_class_species_1 == species]
            species_name = species if species == "Not classified" else "<i>{}</i>".format(species)
            species_full = species_df.iloc[0]["name_classified_species_1"]
            data.append(go.Box(
                go.Box(
                    x=species_df.loc[:, plot_value],
                    text=species_df["name"],
                    marker=dict(color=COLOR_DICT.get(species_full, None)),
                    boxpoints="all",
                    jitter=0.3,
                    pointpos=-1.8,
                    name="{} ({})".format(species_name,species_df["_id"].count()),
                    showlegend=False
                )
            ))
        return {
            "data": data,
            "layout": go.Layout(
                title=plot_value.replace("_", " "),
                margin=go.Margin(
                    l=175,
                    r=50,
                    b=25,
                    t=50
                ),
                xaxis={"showgrid": True}
            )
        }

    @app.callback(
        Output('group-list', 'value'),
        [Input('group-all', 'n_clicks')]
    )
    def all_groups(n_clicks):
        return group_list

    @app.callback(
        Output('species-list', 'value'),
        [Input('species-all', 'n_clicks')]
    )
    def all_species(n_clicks):
        return species_list


    @app.callback(
        Output('sample-list', "value"),
        [Input("summary-plot", "clickData")]
    )
    def update_active_sample(clickData):
        if clickData != None:
            return clickData["points"][0]["text"]

    app.run_server(debug=True)


if __name__ == "__main__":
    main(sys.argv)
