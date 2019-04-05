import dash_html_components as html
import dash_core_components as dcc
import dash_table
import pandas as pd
import numpy as np
import components.import_data as import_data
import components.admin as admin
import dash_bootstrap_components as dbc
import components.global_vars as global_vars

TABLE_PAGESIZE = 100


def format_selected_samples(filtered_df):
    "Returns a formatted string of selected samples"
    return "\n".join([row["name"] for index, row in filtered_df.iterrows()])

def html_div_filter():
    qc_list_options = [
        {"label": "OK", "value": "pass:OK"},
        {"label": "Core Facility", "value": "fail:core facility"},
        {"label": "Supplying Lab", "value": "fail:supplying lab"},
        {"label": "Not checked", "value": "Not checked"}
    ]
    

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
                dbc.ButtonGroup([
                    dbc.Button("Search samples",
                            id="apply-filter-button",
                            color="primary",
                            n_clicks=0),
                    dbc.Button("Generate download link (1000 samples max.)",
                            id="generate-download-button",
                            color="secondary",
                            n_clicks=0)
                ]),
                html.Div(id="tsv-download")
            ]
        ),
        html.Div(
            [
                html.Div([
                    html.H5('Search Results'),
                    html.Div([
                        html.H6([
                            html.Span("0", id="filter-sample-count"),
                            ' samples loaded.'
                        ]),
                        
                    ]),
                    dash_table.DataTable(

                        data=[{}],
                        style_table={
                            'overflowX': 'scroll',
                            'overflowY': 'scroll',
                            'maxHeight': '480'
                        },
                        columns=global_vars.COLUMNS,
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
                        # row_deletable=True,
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
            ]
        )
    ])


def generate_table(tests_df):
    qc_action = "stamps.ssi_stamper.value"
    if qc_action not in tests_df:
        tests_df[qc_action] = np.nan
    else:
        tests_df[qc_action] = tests_df[qc_action].str.split(":", expand=True)[
            1]

    if "reads.R1" not in tests_df:
        tests_df["reads.R1"] = np.nan

    values = {"reads.R1": ""}
    tests_df = tests_df.fillna(value=values)
    no_reads_mask = tests_df["reads.R1"] == ""
    tests_df.loc[no_reads_mask, qc_action] = "core facility (no reads)"
    mask = pd.isnull(tests_df[qc_action])
    tests_df.loc[mask, qc_action] = "not tested"
    slmask = tests_df[qc_action] == "supplying lab"
    tests_df.loc[slmask, qc_action] = "warning: supplying lab"

    user_stamp_col = "stamp.supplying_lab_check.value"
    # Overload user stamp to ssi_stamper
    if user_stamp_col in tests_df.columns:
        user_OK_mask = tests_df[user_stamp_col] == "pass:OK"
        tests_df.loc[user_OK_mask, qc_action] = "*OK"
        user_sl_mask = tests_df[user_stamp_col] == "fail:supplying lab"
        tests_df.loc[user_sl_mask, qc_action] = "*warning: supplying lab"
        user_cf_mask = tests_df[user_stamp_col] == "fail:core facility"
        tests_df.loc[user_cf_mask, qc_action] = "*core facility"

    # Split test columns
    columns = tests_df.columns
    split_columns = [
        "sample_components.ssi_stamper.summary.assemblatron:1x10xsizediff",
        "sample_components.ssi_stamper.summary.whats_my_species:minspecies",
        "sample_components.ssi_stamper.summary.whats_my_species:nosubmitted",
        "sample_components.ssi_stamper.summary.whats_my_species:detectedspeciesmismatch"
    ]
    i = 0
    for column in columns:
        if column in split_columns:
            new = tests_df[column].str.split(":", expand=True)
            loc = tests_df.columns.get_loc(column)
            #tests_df.drop(columns = [column], inplace=True)
            tests_df.insert(loc, column + "_QC", new[0])
            tests_df.insert(loc + 1, column + "_text", new[2])
        i += 1

    if "sample_components.ariba_mlst.summary.mlst_report" in columns:
        first_split = tests_df["sample_components.ariba_mlst.summary.mlst_report"].str.split(
            ",", n=1, expand=True)
        if len(first_split.columns) == 2:
            second_split = first_split[0].str.split(":", n=1, expand=True)
            if len(second_split.columns) == 2:
                keyerrormask = second_split[1] == " 'ariba_mlst/mlst_report_tsv'"
                second_split.loc[keyerrormask, 1] = np.nan
                tests_df["ariba_mlst_type"] = second_split[1]
                tests_df["ariba_mlst_alleles"] = first_split[1]

    test_cols = [col for col in columns if (col.startswith(
        "sample_components.ssi_stamper.summary") and
        not col.startswith("sample_components.ssi_stamper.summary.qcquickie"))]

    def concatenate_failed(row):
        res = []
        for col in test_cols:
            test_name = col.split(":")[-1]
            if type(row[col]) == str:
                fields = row[col].split(":")
                if fields[0] in ["fail", "undefined"]:
                    res.append("Test {}: {}, {}".format(
                        test_name, fields[0], fields[1]))
        row["ssi_stamper_failed_tests"] = ". ".join(res)
        return row

    # Round columns:
    for col in global_vars.ROUND_COLUMNS:
        if col in tests_df.columns:
            tests_df[col] = round(tests_df[col], 3)

    tests_df = tests_df.apply(concatenate_failed, axis="columns")

    COLUMNS = global_vars.COLUMNS

    # Generate conditional formatting:
    style_data_conditional = []
    conditional_columns = [
        col for col in tests_df.columns if col.startswith("QC_")]

    for status, color in ("fail", "#ea6153"), ("undefined", "#f1c40f"):
        style_data_conditional += list(map(lambda x: {"if": {
            "column_id": x, "filter": '{} eq "{}"'.format(x, status)}, "backgroundColor": color}, conditional_columns))

    for status, color in ("core facility", "#ea6153"), ("warning: supplying lab", "#f1c40f"):
        style_data_conditional += [{"if": {
            "column_id": qc_action, "filter": 'QC_action eq "{}"'.format(status)}, "backgroundColor": color}]

    tests_df["_id"] = tests_df["_id"].astype(str)
    tests_df = tests_df.filter([ c["id"] for c in COLUMNS])

    return tests_df
