import dash_html_components as html
import dash_core_components as dcc
import dash_table
import import_data


from components.global_vars import PLOTS, DEFAULT_PLOT


def format_selected_samples(filtered_df):
    "Returns a formatted string of selected samples"
    return "\n".join([row["name"] for index, row in filtered_df.iterrows()])

def html_div_summary():
    plot_values_options = [{"label": plot, "value": plot}
                           for plot, value in PLOTS.items()]
    qc_options = ["OK", "core facility", "supplying lab", "skipped", "Not checked"]
    qc_list_options = [{"label":o, "value":o} for o in qc_options]
    

    run_list = list(import_data.get_run_list())
    run_list_options = [
        {
            "label": "{} ({})".format(run["name"],
                                        len(run["samples"])),
            "value": run["name"]
        } for run in run_list]

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
                                ],
                                id="run-selector",
                                className="row"
                            ),
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
                                            html.Div(
                                                dcc.Dropdown(
                                                    id="group-list",
                                                    multi=True
                                                ),
                                                id="group-div"
                                            )
                                        ],
                                        className="twelve columns",
                                        id="group-form"
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
                                                htmlFor="species-list"
                                            ),
                                            html.Div(
                                                dcc.Dropdown(
                                                    id="species-list",
                                                    multi=True
                                                ),
                                                id="species-div"
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
                                                    "Passed ssi_stamper ",
                                                    html.Small(
                                                        [
                                                            "(",
                                                            html.A(
                                                                "all",
                                                                href="#",
                                                                n_clicks=0,
                                                                id="qc-all"
                                                            ),
                                                            ")"
                                                        ]
                                                    )
                                                ],
                                                htmlFor="qc-list"
                                            ),
                                            html.Div(
                                                dcc.Dropdown(
                                                    id="qc-list",
                                                    multi=True,
                                                    options=qc_list_options,
                                                    value=qc_options
                                                ),
                                                id="qc-div"
                                            )
                                        ],
                                        className="twelve columns"
                                    )
                                ],
                                className="row"
                            )
                        ],
                        className="twelve columns"
                    )
                ], className="row"),
                html.Div([
                    html.Div(
                        [
                            html.Button(
                                "Apply Filter",
                                id="apply-filter-button",
                                n_clicks=0,
                                className="button-primary u-full-width"
                            )
                        ],
                        className="twelve columns"
                    ),
                    html.Div(
                        [
                            html.Label(
                                [
                                    "Selected Samples ():"
                                ],
                                htmlFor="plot-list"),
                            dcc.Textarea(
                                className="u-full-width",
                                style={"resize": "none",
                                       "height": "300px"},
                                readOnly=True,
                                value=[""],
                                id="selected-samples-list"
                            ),
                            html.Div("",
                                     style={"display": "none"},
                                     id="selected-samples-ids")
                        ],
                        style={"display": "none"},
                        #className="four columns",
                        id="selected-samples"
                    ),
                    
                ],
                className="row mt-1"
            ),
            html.Div([
                    html.H6('No samples loaded. Click "Apply Filter" to load samples.'),
                    html.Div([
                        html.P(
                            'To filter on a string type eq, space and exact text in double quotes: eq "FBI"'),
                        html.P(
                            'To filter on a number type eq, < or >, space and num(<number here>): > num(500)')
                    ]),
                    html.Div(dash_table.DataTable(id="datatable-ssi_stamper", data=[{}]), style={"display": "none"})
                ],
                     id="ssi_stamper-report", className="bigtable"),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Plot value",
                                       htmlFor="plot-list"),
                            dcc.Dropdown(
                                id="plot-list",
                                options=plot_values_options,
                                value=DEFAULT_PLOT
                            )
                        ],
                        className="twelve columns"
                    )
                ],
                className="row"
            ),
            dcc.Graph(id="summary-plot"),
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            html.Div("",
                                                     style={"display": "none"},
                                                     id="lasso-sample-ids")
                                        ],
                                        className="twelve columns",
                                        id="lasso-div"
                                    )
                                ],
                                className="row"
                            )
                        ],
                        className="twelve columns"
                    )
                ]
                , className="row")
        ], className="border-box"
    )
