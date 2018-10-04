import dash_html_components as html
import dash_core_components as dcc

from components.global_vars import PLOTS, DEFAULT_PLOT


def format_selected_samples(filtered_df):
    "Returns a formatted string of selected samples"
    return "\n".join([row["name"] for index, row in filtered_df.iterrows()])

def html_div_summary():
    plot_values_options = [{"label": plot, "value": plot}
                           for plot, value in PLOTS.items()]
    return html.Div(
        [
            html.H5("Summary", className="box-title"),
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
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
                            )
                        ],
                        className="eight columns"
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
                        className="four columns",
                        id="selected-samples"
                    ),
                    
                ],
                className="row"
            ),
            dcc.Graph(id="summary-plot"),
            html.Div(id="testomatic-report", className="bigtable"),
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
