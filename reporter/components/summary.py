import dash_html_components as html
import dash_core_components as dcc
import dash_table
import components.import_data as import_data
import components.admin as admin


filter_notice_div = html.Div([
    html.P(
        'To filter on a string type eq, space and exact text in double quotes: eq "FBI"'),
    html.P(
        'To filter on a number type eq, < or >, space and num(<number here>): > num(500)'),
    html.P('* on QC: manually curated')
])


def format_selected_samples(filtered_df):
    "Returns a formatted string of selected samples"
    return "\n".join([row["name"] for index, row in filtered_df.iterrows()])

def html_div_summary():
    qc_options = ["OK", "core facility", "supplying lab", "skipped", "Not checked"]
    qc_list_options = [{"label":o, "value":o} for o in qc_options]
    

    return html.Div([
        html.Div(
            [
                html.H5("Filter", className="box-title"),
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
                                                        options=[],
                                                        placeholder="Sequencing run"
                                                    ),
                                                    id="run-list-div"
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
                                                        "Supplying lab"
                                                    ],
                                                    htmlFor="group-list"
                                                ),
                                                html.Div(
                                                    dcc.Dropdown(
                                                        id="group-list",
                                                        multi=True,
                                                        value=[]
                                                    ),
                                                    id="group-div"
                                                )
                                            ],
                                            className="six columns",
                                            id="group-form"
                                        ),
                                        html.Div(
                                            [
                                                html.Label(
                                                    [
                                                        "Species: "
                                                    ],
                                                    htmlFor="species-list",
                                                    style={"display": "inline-block"}
                                                ),
                                                dcc.RadioItems(
                                                    options=[
                                                        {"label": "Provided",
                                                        "value": "provided"},
                                                        {"label": "Detected",
                                                        "value": "detected"},
                                                    ],
                                                    value="provided",
                                                    labelStyle={
                                                        'display': 'inline-block'},
                                                    id="form-species-source",
                                                    style={
                                                        "display": "inline-block"}

                                                ),
                                                html.Div(
                                                    dcc.Dropdown(
                                                        id="species-list",
                                                        multi=True,
                                                        value=[]
                                                    ),
                                                    id="species-div"
                                                )
                                            ],
                                            className="six columns"
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
                                                        "Passed ssi_stamper"
                                                    ],
                                                    htmlFor="qc-list"
                                                ),
                                                html.Div(
                                                    dcc.Dropdown(
                                                        id="qc-list",
                                                        multi=True,
                                                        options=qc_list_options,
                                                        placeholder="All values selected",
                                                        value=[]
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
                    ], className="row"
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    html.Button(
                                        "Apply Filter",
                                        id="apply-filter-button",
                                        n_clicks=0,
                                        n_clicks_timestamp=0,
                                        className="button-primary u-full-width"
                                    ), id="applybutton-div"
                                )
                            ],
                            className="twelve columns"
                        ),
                    ],
                    className="row mt-1"
                ),
            ], className="border-box"
        ),
        html.Div(
            [
                html.H5("Summary", className="box-title"),
                html.Div([
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
                    html.Div(dash_table.DataTable(id="datatable-ssi_stamper", data=[{}]), style={"display": "none"})
                ], id="ssi_stamper-report", className="bigtable"),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Label("Plot species",
                                        htmlFor="plot-species"),
                                dcc.RadioItems(
                                    options=[
                                        {"label": "Provided",
                                        "value": "provided"},
                                        {"label": "Detected",
                                        "value": "detected"},
                                    ],
                                    value="provided",
                                    labelStyle={
                                        'display': 'inline-block'},
                                    id="plot-species-source"
                                ),
                                html.Div(dcc.Dropdown(
                                    id="plot-species"
                                ),id="plot-species-div")
                                
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
                    , className="row"),
                html.Div(admin.selected_samples_div(),className="row")
            ], className="border-box"
        )
    ])
