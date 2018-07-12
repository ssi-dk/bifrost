# -*- coding: utf-8 -*-
import os
import sys

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
from components.summary import html_div_summary, format_selected_samples
from components.sample_report import children_sample_list_report
from components.images import list_of_images, static_image_route, get_species_color, COLOR_DICT, image_directory
import components.global_vars


def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

# Globals

PAGESIZE = 50

def filter_dataframe(dataframe, species_list, group_list, run_name, page_n=None):
    if species_list is None: species_list = []
    if group_list is None: group_list = []
    filtered = dataframe[
        dataframe.short_class_species_1.isin(species_list) &
        dataframe.supplying_lab.isin(group_list)
    ]
    if len(run_name) != 0:
        filtered = filtered[filtered.run_name == run_name]
    if page_n is None:
        return filtered
    else:
        return filtered.iloc[page_n*PAGESIZE:(page_n+1)*PAGESIZE]

def paginate_df(dataframe, page_n):
    return dataframe.iloc[page_n*PAGESIZE:(page_n+1)*PAGESIZE]

def html_table_run_information(run_data):
    return html_table([
        ["Run Name", run_data["run_name"]],
        ["Run Date", "06 MAY 2018"],
        ["Placeholder", "Lorem ipsum"]
    ])

def main(argv):
    dataframe = import_data.import_data()
    

    # Change this to provided species ?
    dataframe = dataframe.sort_values("short_class_species_1")
    app = dash.Dash()
    app.config['suppress_callback_exceptions'] = True

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
        dash_scroll_up.DashScrollUp(
            id='input',
            label='UP',
            className="button button-primary no-print"
        ),
        html.Div(dt.DataTable(rows=[{}], editable=False), style={'display': 'none'}),
        html.H1("QC REPORT"),
        html.H2("", id="run-name"),
        dcc.Location(id="url", refresh=False),
        html.Div(html_table_run_information({'run_name': ""}), id="run-table"),
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
                                    "QCQuickie report",
                                    id="update-qcquickie",
                                    n_clicks_timestamp=0,
                                    className="button-primary u-full-width"
                                )
                            ],
                            className="four columns"
                        ),
                        html.Div(
                            [
                                html.Button(
                                    "Assembly report",
                                    id="update-assembly",
                                    n_clicks_timestamp=0,
                                    className="button-primary u-full-width"
                                )
                            ],
                            className="four columns"
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
                            className="four columns"
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
        Output('run-name', 'children'),
        [Input('url', 'pathname')]
    )
    def update_run_name(pathname):
        if pathname is None:
            pathname = "/"
        path = pathname.split('/')
        if path[1] in dataframe.run_name.values:
            return path[1]
        elif path[1] == "":
            return []
        else:
            return "Not found"

    @app.callback(
        Output('run-selector', 'children'),
        [Input('run-name', 'children')]
    )
    def update_run_list(run_name):
        if len(run_name) == 0:
            run_list = dataframe.sort_values("run_name").run_name.unique()
            run_list_options = [{"label": run, "value": run}
                                for run in run_list]
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
        Output('report-count', 'children'),
        [Input('summary-plot', 'selectedData'),
         Input(component_id="species-list", component_property="value"),
         Input(component_id="group-list", component_property="value"),
         Input(component_id="run-name", component_property="children")]
    )
    def display_selected_data(selected_data, species_list, group_list, run_name):
        if selected_data is not None and len(selected_data["points"]):
            return len([sample['text']
                      for sample in selected_data["points"]])
        else:
            filtered_df = filter_dataframe(
                dataframe, species_list, group_list, run_name)
            return filtered_df['name'].count()

    @app.callback(
        Output('lasso-div', 'children'),
        [Input('summary-plot', 'selectedData'),
         Input('report-count', 'children')]
    )
    def display_selected_data(selected_data, ignore_this):
        # ignore_this is there so the function is called when the sample list is updated.
        if selected_data is not None and len(selected_data["points"]):
            points = [sample['text']
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
                )
            ]
        else:
            return ""

    @app.callback(
        Output('group-div', 'children'),
        [Input('run-name', 'children')]
    )
    def update_group_list(run_name):
        if len(run_name) == 0:
            group_list = dataframe.supplying_lab.value_counts()
        else:
            group_list = dataframe[dataframe.run_name == run_name].supplying_lab.value_counts()
        group_list_options = [{"label": "{} ({})".format(group, count), "value": group}
                              for group, count in group_list.items()]
        group_options = [group for group, count in group_list.items()]
        return dcc.Dropdown(
            id="group-list",
            options=group_list_options,
            multi=True,
            value=group_options
        )

    @app.callback(
        Output('species-div', 'children'),
        [Input('run-name', 'children')]
    )
    def update_group_list(run_name):
        if len(run_name) == 0:
            species_list = dataframe.short_class_species_1.value_counts()
        else:
            species_list = dataframe[dataframe.run_name ==
                                     run_name].short_class_species_1.value_counts()
        species_list_options = [{"label": "{} ({})".format(species, count), "value": species}
                                for species, count in species_list.items()]
        species_options = [species for species, count in species_list.items()]
        return dcc.Dropdown(
            id="species-list",
            options=species_list_options,
            multi=True,
            value=species_options
        )


    @app.callback(
        Output('run-table', 'children'),
        [Input('run-name', 'children')]
    )
    def update_run_name(run_name):
        if run_name == "Not found":
            return html_table_run_information({'run_name': "Run not found!"})
        elif run_name == None:
            return html_table_run_information({'run_name': "No run selected"})
        else:
            return html_table_run_information({'run_name': run_name})

    @app.callback(
        Output(component_id="page-n",
               component_property="children"),
        [Input(component_id="prevpage", component_property="n_clicks_timestamp"),
         Input(component_id="nextpage", component_property="n_clicks_timestamp")],
        [State(component_id="page-n", component_property="children"),
         State(component_id="species-list", component_property="value"),
         State(component_id="group-list", component_property="value"),
         State(component_id="run-name", component_property="children")]
    )
    def next_page(prev_ts, next_ts, page_n, species_list, group_list, run_name):
        page_n = int(page_n)
        max_page = filter_dataframe(dataframe, species_list, group_list, run_name)["_id"].count() // PAGESIZE
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
         [State('summary-plot', 'selectedData'),
         State("selected-samples", "value")]
         )
    def sample_report(page_n, data_content, lasso_selected, prefilter_samples):
        if not (data_content == "qcquickie" or data_content == "assembly"):
            return []
        page_n = int(page_n)
        if lasso_selected is not None and len(lasso_selected["points"]):
            samples = [sample['text']
                       for sample in lasso_selected["points"]]  # lasso first
        else:
            samples = prefilter_samples[0].split("\n")
        filtered = dataframe[dataframe.name.isin(samples)]
        page = paginate_df(filtered, page_n)
        max_page = len(filtered) // PAGESIZE
        page_species = page.qcquickie_name_classified_species_1.unique().tolist()
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
        [Input(component_id="page-n", component_property="children")],
        [State(component_id="species-list", component_property="value"),
         State(component_id="group-list", component_property="value"),
         State(component_id="run-name", component_property="children")]
    )
    def update_nextpage(page_n, species_list, group_list, run_name):
        page_n = int(page_n)
        max_page = len(filter_dataframe(
            dataframe, species_list, group_list, run_name)) // PAGESIZE
        if page_n == max_page:
            return True
        else:
            return False

    @app.callback(
        Output(component_id="current-report", component_property="children"),
        [Input(component_id="update-qcquickie", component_property="n_clicks_timestamp"),
         Input(component_id="update-assembly", component_property="n_clicks_timestamp"),
         Input(component_id="update-table", component_property="n_clicks_timestamp")],
        [State('summary-plot', 'selectedData'),
         State("selected-samples", "value")]
    )
    def update_report(n_qcquickie_ts, n_assembly_ts, n_table_ts, lasso_selected, prefilter_samples):
        if max(n_qcquickie_ts, n_assembly_ts) > n_table_ts:  # samples was clicked
            if n_qcquickie_ts > n_assembly_ts:
                title = "QCQuickie Report"
                content = "qcquickie"
            else:
                title = "Assembly Report"
                content = "assembly"
            return [
                html.H3(title),
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
                
                html.Div(id="sample-report", **{"data-content": content}),
            ]
        elif n_table_ts > n_assembly_ts:  # table was clicked
            if lasso_selected is not None and len(lasso_selected["points"]):
                samples = [sample['text'] for sample in lasso_selected["points"]] # lasso first
            else:
                samples = prefilter_samples[0].split("\n")
            
            return [
                html.H3("Table Report"),

                dt.DataTable(
                    rows=dataframe[dataframe.name.isin(samples)].to_dict("records"),

                    columns=global_vars.columns, # sets the order
                    column_widths=[150]*len(global_vars.columns),
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
         Input(component_id="group-list", component_property="value"),
         Input(component_id="run-name", component_property="children")]
    )
    def update_selected_samples(species_list, group_list, run_name):
        filtered_df = filter_dataframe(dataframe, species_list, group_list, run_name)
        return [format_selected_samples(filtered_df)]

    @app.callback(
        Output(component_id="sample-count", component_property="children"),
        [Input(component_id="species-list", component_property="value"),
         Input(component_id="group-list", component_property="value"),
         Input(component_id="run-name", component_property="children")]
    )
    def update_sample_count(species_list, group_list, run_name):
        filtered_df = filter_dataframe(dataframe, species_list, group_list, run_name)
        return filtered_df['name'].count()

    @app.callback(
        Output(component_id="summary-plot", component_property="figure"),
        [Input(component_id="species-list", component_property="value"),
         Input(component_id="group-list", component_property="value"),
         Input(component_id="run-name", component_property="children"),
         Input(component_id="plot-list", component_property="value")]
    )
    def update_coverage_figure(species_list, group_list, run_name, plot_value):
        data = []
        filtered_df = filter_dataframe(dataframe, species_list, group_list, run_name)
        if species_list is None: species_list = []
        for species in species_list:
            species_df = filtered_df[filtered_df.short_class_species_1 == species]
            if species == "Not classified":
                species_name = species
            else:
                species_name = "<i>{}</i>".format(species)
            data.append(go.Box(
                go.Box(
                    x=species_df.loc[:, plot_value],
                    text=species_df["name"],
                    marker=dict(
                        color=COLOR_DICT.get(species, None),
                        size=4
                    ),
                    boxpoints="all",
                    jitter=0.5,
                    pointpos=-1.8,
                    name="{} ({})".format(species_name,species_df["_id"].count()),
                    showlegend=False
                )
            ))
        return {
            "data": data,
            "layout": go.Layout(
                hovermode="closest",
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
        [Input('group-all', 'n_clicks')],
        [State("run-name", "children")]
    )
    def all_groups(n_clicks, run_name):
        return filter_dataframe(dataframe, [], [], run_name).supplying_lab.unique()

    @app.callback(
        Output('species-list', 'value'),
        [Input('species-all', 'n_clicks')],
        [State("run-name", "children")]
    )
    def all_species(n_clicks, run_name):
        return filter_dataframe(dataframe, [], [], run_name).short_class_species_1.unique()


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
