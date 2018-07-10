import dash_html_components as html
import dash_core_components as dcc

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


def format_selected_samples(filtered_df):
    "Returns a formatted string of selected samples"
    return "\n".join([row["name"] for index, row in filtered_df.iterrows()])

def html_div_summary():
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
                                                    html.Span(
                                                        id="sample-count"),
                                                    "):"
                                                ],
                                                htmlFor="plot-list"),
                                            dcc.Textarea(
                                                id="selected-samples",
                                                className="u-full-width",
                                                style={"resize": "none",
                                                       "height": "300px"},
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
                    ),
                    
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
