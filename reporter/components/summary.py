import dash_html_components as html
import dash_core_components as dcc

PLOT_VALUES = [
    "qcquickie_bin_length_at_1x",
    "qcquickie_bin_length_at_10x",
    "qcquickie_bin_length_at_25x",
    "qcquickie_bin_length_1x_25x_diff",
    "qcquickie_bin_contigs_at_1x",
    "qcquickie_bin_contigs_at_10x",
    "qcquickie_bin_contigs_at_25x",
    "qcquickie_N50",
    "qcquickie_N75",
    "assembly_bin_length_at_10x",
    "assembly_bin_length_at_1x",
    "assembly_bin_length_at_25x",
    "assembly_bin_length_1x_25x_diff",
    "assembly_bin_contigs_at_1x",
    "assembly_bin_contigs_at_10x",
    "assembly_bin_contigs_at_25x",
    "assembly_N50",
    "assembly_N75"
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
                        className="four columns",
                        id="selected-samples"
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
