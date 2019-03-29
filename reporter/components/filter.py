import dash_html_components as html
import dash_core_components as dcc
import dash_table
import components.import_data as import_data
import components.admin as admin
import dash_bootstrap_components as dbc
import components.global_vars as global_vars

TABLE_PAGESIZE = 100


def simplify_name(name):
    return name.replace(":", "_").replace(".", "_").replace("=", "_")

COLUMNS = []
for column in global_vars.COLUMNS:
    column["id"] = simplify_name(column["id"])
    COLUMNS.append(column)


def format_selected_samples(filtered_df):
    "Returns a formatted string of selected samples"
    return "\n".join([row["name"] for index, row in filtered_df.iterrows()])

def html_div_summary():
    qc_options = ["OK", "core facility", "supplying lab", "skipped", "Not checked"]
    qc_list_options = [{"label":o, "value":o} for o in qc_options]
    

    return html.Div([
        html.Div(
            [
                dcc.Store(id="param-store", storage_type="memory", data={}),
                html.H1("Filter", className="mt-3"),

                dbc.Row([
                    dbc.Col([
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.FormGroup(
                                        [
                                            dbc.Label(
                                                "Run", html_for="run-list"),
                                            dcc.Dropdown(
                                                id="run-list",
                                                multi=True,
                                                value=[],
                                                placeholder="All runs selected",
                                            ),
                                        ]
                                    ),
                                    width=6,
                                ),
                                dbc.Col(
                                    dbc.FormGroup(
                                        [
                                            dbc.Label(
                                                "Group", html_for="group-list"),
                                            dcc.Dropdown(
                                                id="group-list",
                                                multi=True,
                                                value=[],
                                                placeholder="All groups selected",
                                            ),
                                        ]
                                    ),
                                    width=6,
                                ),
                            ],
                            form=True,
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.FormGroup(
                                        [
                                            html.Label(
                                                [
                                                    "Species: "
                                                ],
                                                htmlFor="species-list",
                                                style={
                                                    "display": "inline-block"}
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
                                                    value=[],
                                                    placeholder="All species selected",
                                                ),
                                                id="species-div"
                                            )
                                        ]
                                    ),
                                    width=6,
                                ),
                                dbc.Col(
                                    dbc.FormGroup(
                                        [
                                            dbc.Label(
                                                "SSI QC", html_for="qc-list"),
                                            dcc.Dropdown(
                                                id="qc-list",
                                                multi=True,
                                                options=qc_list_options,
                                                placeholder="All values selected",
                                                value=[]
                                            ),
                                        ]
                                    ),
                                    width=6,
                                ),
                            ],
                            form=True,
                        ),
                    ], width=9),
                    dbc.Col([
                        dbc.Row(dbc.Col([
                            dbc.Label("Sample names", html_for="samples-form"),
                            dcc.Textarea(
                                id="samples-form",
                                placeholder="one sample per line",
                                value="",
                                rows=6
                            )
                        ], width=12))
                    ], width=3),
                ]),
                dbc.Button("Search samples",
                           id="apply-filter-button",
                           color="primary",
                           n_clicks=0,
                           className="mr-1"),
                dbc.Button("Generate download link",
                           color="secondary",
                           className="mr-1"),
            ]
        ),
        html.Div(
            [
                html.Div([
                    html.H6('No samples loaded. Click "Apply Filter" to load samples.'),
                    html.Div([
                        html.Div([
                            "Download Table ",
                            html.A("(tsv, US format)",
                                id="ustsv-download",
                                download='report.tsv'),
                            " - ",
                            html.A("(csv, EUR Excel format)",
                                id="eurtsv-download",
                                download='report.csv')
                        ], className="six columns"),
                        
                    ], className="row"),
                    dash_table.DataTable(

                        data=[{}],
                        style_table={
                            'overflowX': 'scroll',
                            'overflowY': 'scroll',
                            'maxHeight': '480'
                        },
                        columns=COLUMNS,
                        style_cell={
                            'width': '200px',
                            'padding': '0 15px'
                        },
                        style_cell_conditional=[
                            {
                                "if": {"column_id": "ssi_stamper_failed_tests"},
                                "textAlign": "left"
                            }
                        ],
                        n_fixed_rows=1,
                        row_deletable=True,
                        # filtering=True,  # Front end filtering
                        # sorting=True,
                        selected_rows=[],
                        # style_data_conditional=style_data_conditional,
                        pagination_settings={
                            'current_page': 0,
                            'page_size': TABLE_PAGESIZE
                        },
                        pagination_mode='be',
                        id="datatable-ssi_stamper"
                    )
                ], id="ssi_stamper-report", className="bigtable"),
                # html.Div(
                #     [
                #         html.Div(
                #             [
                #                 html.Label("Plot species",
                #                         htmlFor="plot-species"),
                #                 dcc.RadioItems(
                #                     options=[
                #                         {"label": "Provided",
                #                         "value": "provided"},
                #                         {"label": "Detected",
                #                         "value": "detected"},
                #                     ],
                #                     value="provided",
                #                     labelStyle={
                #                         'display': 'inline-block'},
                #                     id="plot-species-source"
                #                 ),
                #                 html.Div(dcc.Dropdown(
                #                     id="plot-species"
                #                 ),id="plot-species-div")
                                
                #             ],
                #             className="twelve columns"
                #         )
                #     ],
                #     className="row"
                # ),
                # dcc.Graph(id="summary-plot"),
                # html.Div(
                #     [
                #         html.Div(
                #             [
                #                 html.Div(
                #                     [
                #                         html.Div(
                #                             [
                #                                 html.Div("",
                #                                         style={"display": "none"},
                #                                         id="lasso-sample-ids")
                #                             ],
                #                             className="twelve columns",
                #                             id="lasso-div"
                #                         )
                #                     ],
                #                     className="row"
                #                 )
                #             ],
                #             className="twelve columns"
                #         )
                #     ]
                #     , className="row"),
            ], className="border-box"
        )
    ])
