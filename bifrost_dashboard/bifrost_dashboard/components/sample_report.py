import dash_html_components as html
import dash_core_components as dcc
import dash_table as dt
from bifrost_dashboard.components.images import list_of_images
from bifrost_dashboard.components.table import html_table, html_td_percentage
import bifrost_dashboard.components.import_data as import_data
import bifrost_dashboard.components.global_vars as global_vars
import bifrost_dashboard.components.admin as admin
import dash_bootstrap_components as dbc
import pandas as pd
import numpy as np
import math
import json

SAMPLE_PAGESIZE = 10

def get(dict_o, field, default=None):
    r = dict_o.get(field, default)
    
    if not isinstance(r, (list, np.ndarray)) and pd.isna(r):
        return default
    else:
        return r

def sample_report(data):

    sample_n = len(data)

    return [
        html.Span("0", style={"display": "none"}, id="page-n"),
        html.Span(str(sample_n // SAMPLE_PAGESIZE),
                  style={"display": "none"}, id="max-page"),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Button(
                            "Previous page",
                            id="prevpage",
                            className="btn btn-outline-secondary",
                            n_clicks_timestamp=0),
                    ],
                    width="auto"
                ),
                dbc.Col(
                    [
                        html.Button(
                            "Next page",
                            id="nextpage",
                            className="btn btn-outline-secondary",
                            n_clicks_timestamp=0),
                    ],
                    width="auto"
                ),
            ], justify="between"
        ),

        html.Div(id="sample-report"),

        html.Div(dbc.Row(
            [
                dbc.Col(
                    [
                        html.Button(
                            "Previous page",
                            id="prevpage2",
                            className="btn btn-outline-secondary",
                            n_clicks_timestamp=0),
                    ],
                    width="auto"
                ),
                dbc.Col(
                    [
                        html.Button(
                            "Next page",
                            id="nextpage2",
                            className="btn btn-outline-secondary",
                            n_clicks_timestamp=0),
                    ],
                    width="auto"
                ),
            ], justify="between"
        ), className="mb-4"),
    ]

def check_test(test_name, sample):
    test_path = "properties.stamper.summary." + test_name + ".status"
    res = get(sample, test_path, "")
    if res == "":
        return "" # show nothing
        #return "test-missing"
    if res.startswith("pass"):
        return "test-pass"
    elif res.startswith("fail"):
        return "test-fail"
    else:
        return "test-warning"

def get_species_img(sample_data):
    genus = str(get(sample_data, "properties.species_detection.summary.name_classified_species_1")).split()[
        0].lower()
    if "{}.svg".format(genus) in list_of_images:
        img = html.Img(
            src="/assets/img/" + genus + ".svg",
            className="svg_bact"
        )
    else:
        img = None
    return img

def generate_sample_report(sample, n_sample):
    img = get_species_img(sample)
    img_div = html.Div(
            [
                img
            ],
            className="d-inline")
    return (
        dbc.Card(
            [
                dbc.CardHeader([
                    img_div,
                    html.Div([
                        html.A(id="sample-" + str(sample["name"])),
                        html.H6(
                            sample["name"],
                            className="d-inline font-weight-bold text-primary mx-2"),
                        html.I(get(sample, "properties.species_detection.summary.species"))
                    ], className="d-inline-block"),
                    
                ]),
                dbc.CardBody([
                    html_sample_tables(sample),
                    admin.sample_radio_feedback(sample, n_sample)
                ])
                
            ], className="shadow mb-4"
        )
    )


def html_organisms_table(sample_data, **kwargs):
    percentages = [
        get(sample_data, "properties.species_detection.summary.percent_classified_species_1", math.nan),
        get(sample_data, "properties.species_detection.summary.percent_classified_species_2", math.nan),
        get(sample_data, "properties.species_detection.summary.percent_unclassified", math.nan)
    ]

    color_0 = "#b3ccc1"
    color_1 = "#b3ccc1" #Green

    color_2 = "#f3bbd3"  # Default

#   color_u = COLOR_DICT.get("", "#fee1cd")  # Default
    color_u = "#fee1cd"  # Default

    return html.Div([
        html.Table([
            html.Tr([
                html.Td(
                    [html.I(get(sample_data, "properties.species_detection.summary.name_classified_species_1", "No data")), " + Unclassified"], className="cell"),
                html_td_percentage(percentages[0] + percentages[2], color_0)
            ], className=check_test("test__species_detection__main_species_level", sample_data) + " trow"),
            html.Tr([
                html.Td(
                    html.I(get(sample_data, "properties.species_detection.summary.name_classified_species_1", "No data")), className="cell"),
                html_td_percentage(percentages[0], color_1)
            ], className="trow"),
            html.Tr([
                html.Td(
                    html.I(get(sample_data, "properties.species_detection.summary.name_classified_species_2", "No data")), className="cell"),
                html_td_percentage(percentages[1], color_2)
            ], className="trow"),
            html.Tr([
                html.Td("Unclassified", className="cell"),
                html_td_percentage(percentages[2], color_u)
            ], className=check_test("test__species_detection__unclassified_level", sample_data) + " trow")
        ], className="bifrost-table")
    ], **kwargs)

def html_test_tables(sample_data):
    qc_to_check = ["stamper"]
    rows = []
    for key, value in sample_data.items():
        if key.startswith("properties.stamper.summary.whats_my_species") \
                or key.startswith("properties.stamper.summary.assemblatron"):
            if pd.isnull(value):
                value = "nan"
            name_v = key.split(":")
            values = value.split(':')
            if len(values) == 3:
                if values[0] != "pass":
                    rows.append([name_v[1].capitalize(), values[1]])
            if (key.endswith(".action")):
                rows.append(["QC Action", value])

    qc_rows = []
    for qc_cat in qc_to_check:
        qc_key = "properties.{}.summary.stamp".format(qc_cat)
        qc_key_val = "{}.value".format(qc_key)
        qc_val = get(sample_data, qc_key_val)
        if qc_val is not None:
            if str(qc_val).startswith("OK"):
                qc_class = "test-pass"
            elif str(qc_val) in ["core facility", "supplying lab", "other"]:
                qc_class = "test-fail"
            else:
                qc_class = ""
            source = get(sample_data, "{}.display_name".format(qc_key))
            reason = get(sample_data, "{}.reason".format(qc_key), "")
            qc_label = "{} [{}, {}]".format(qc_cat, source, reason)
            qc_rows.append({
                "list": [qc_label, qc_val],
                "className": qc_class
            })

    if len(rows):
        test_table = html_table(rows)
    else:
        test_table = html.P("No failed tests.")
    
    if len(qc_rows):
        qc_table = html_table(qc_rows)
    else:
        qc_table = html.P("No info.")

    return dbc.Row([
        dbc.Col([
            html.H6("QC", className="table-header"),
            qc_table
        ]),
        dbc.Col([
            html.H6("Failed QC tests", className="table-header"),
            test_table
        ]),
    ])


def html_sample_tables(sample_data, **kwargs):
    """Generate the tables for each sample containing submitter information,
       detected organisms etc. """

    if get(sample_data, "properties.sample_info.summary.sample_name"):
        emails_l = get(sample_data, "properties.sample_info.summary.emails")
        if emails_l:
            n_emails = len(emails_l.split(";"))
            if (n_emails > 1):
                emails = ", ".join(emails_l.split(";")[:2])
                if (n_emails > 2):
                    emails += ", ..."
            else:
                emails = emails_l
        else:
            emails = ""
        sample_sheet_table = html_table([
                    ["Supplied name", get(sample_data, "sample_sheet.sample_name", "")],
                    ["User Comments", get(sample_data, "sample_sheet.Comments", "")],
                    ["Supplying lab", get(sample_data, "sample_sheet.group", "")],
                    ["Submitter emails", emails],
                    {
                        "list": ["Provided species", html.I(
                            get(sample_data, "sample_sheet.provided_species", ""))],
                        "className": check_test("test__sample__species_provided_is_detected", sample_data)
                    },
                    {
                        "list": ["Read file", 
                                 str(get(sample_data, "properties.datafiles.summary.paired_reads", [""])[0]).split("/")[-1]],
                        "className": check_test("test__sample__has_reads_files", sample_data)
                    }
                ])
    else:
        sample_sheet_table = html_table([])

    title = "Assemblatron Results"
    table = html.Div([
        html_table([
            {
                "list": [
                    "Number of filtered reads",
                    "{:,.0f}".format(
                        get(sample_data, "properties.denovo_assembly.summary.filtered_reads_num", math.nan))
                ],
                "className": check_test("test__denovo_assembly__minimum_read_number", sample_data)
            },
            [
                "Number of contigs (1x cov.)",
                "{:,.0f}".format(
                    get(sample_data, "properties.denovo_assembly.summary.bin_contigs_at_1x", math.nan))
            ],
            [
                "Number of contigs (10x cov.)",
                "{:,.0f}".format(
                    get(sample_data, "properties.denovo_assembly.summary.bin_contigs_at_10x", math.nan))
            ],
            [
                "N50",
                "{:,}".format(
                    get(sample_data, "properties.denovo_assembly.summary.N50", math.nan))
            ],
            {
                "list": [
                    "Average coverage (1x)",
                    "{:,.2f}".format(
                        get(sample_data, "properties.denovo_assembly.summary.bin_coverage_at_1x", math.nan))
                ],
                "className": check_test("test__denovo_assembly__genome_average_coverage", sample_data)
            },
            {
                "list": [
                    "Genome size at 1x depth",
                    "{:,.0f}".format(
                        get(sample_data, "properties.denovo_assembly.summary.bin_length_at_1x", math.nan))
                ],
                "className": check_test("test__denovo_assembly__genome_size_at_1x", sample_data)
            },
            {
                "list": [
                    "Genome size at 10x depth",
                    "{:,.0f}".format(
                        get(sample_data, "properties.denovo_assembly.summary.bin_length_at_10x", math.nan))
                ],
                "className": check_test("test__denovo_assembly__genome_size_at_10x", sample_data)
            },
            {
                "list": [
                    "Genome size 1x - 10x diff",
                    "{:,.0f}".format(
                        get(sample_data, "properties.denovo_assembly.summary.bin_length_at_1x", math.nan)
                        - get(sample_data, "properties.denovo_assembly.summary.bin_length_at_10x", math.nan)
                    )
                ],
                "className": check_test("test__denovo_assembly__genome_size_difference_1x_10x", sample_data)
            },
            [
                "Genome size at 25x depth",
                "{:,.0f}".format(
                    get(sample_data, "properties.denovo_assembly.summary.bin_length_at_25x", math.nan))
            ],
            [
                "Ambiguous sites",
                "{:,.0f}".format(
                    get(sample_data, "properties.denovo_assembly.summary.snp_filter_10x_10%", math.nan))
            ] 
        ])
    ])
    expected_results = global_vars.expected_results
    any_results = False
    results = []
    for entry in expected_results:
        entry_key = "report.{}.data".format(entry)
        if entry_key in sample_data:
            any_results = True
            if len(get(sample_data, entry_key, [])) == 0:
                datatable = "No results"
            else:
                data = get(sample_data, "report.{}.data".format(entry), [])
                columns = get(sample_data, "report.{}.columns".format(entry), [])
                header = [c["name"] for c in columns]
                rows = []
                for rowd in data:

                    row = [rowd[col["id"]] for col in columns]
                    rows.append(row)

                datatable = html_table(rows, header)
                # dt.DataTable(
                #     style_table={
                #         'overflowX': 'scroll',
                #         'overflowY': 'scroll',
                #         'maxHeight': '480'
                #     },

                #     columns=get(
                #         sample_data, "report.{}.columns".format(entry), []),
                #     data=get(sample_data, "report.{}.data".format(entry), []),
                #     page_action='none'
                # )
            results.append(html.Div([
                html.H6(get(sample_data, "report.{}.title".format(entry)),
                        className="table-header"),
                html.Div(get(sample_data, "report.{}.info".format(entry))),
                html.Div(
                    datatable,
                    className="grey-border")
            ], className="col-6"))
        else:
            results.append(html.Div([
                html.H6("{} not run".format(entry.capitalize()),
                        className="table-header")
            ], className="col-6"))

    # mlst_db = sample_data.get("ariba_mlst.mlst_db", "")

    if any_results:
        res_div = html.Details([
            html.Summary(
                "ResFinder/PlasmidFinder/VirulenceFinder/MLST (click to show)"),
            html.Div([
                html.Div(results, className="row"),
            ])
        ])
    else:
        res_div = html.Div(
            html.P("Resfinder, plasmidfinder, MLST and virulencefinder were not run."))


    mlst_type = "ND"
    strain_list = get(sample_data, "properties.mlst.summary.strain")
    if isinstance(strain_list, (list, np.ndarray)):
        mlst_type = ", ".join(list(map(str, strain_list)))


    return html.Div([
        dbc.Row([
            dbc.Col([
                html.H6("Sample Sheet", className="table-header"),
                sample_sheet_table,
                html.H6("Detected Organisms", className="table-header"),
                html_organisms_table(sample_data)
            ]),
            dbc.Col([
                html.H6(title, className="table-header"),
                table,
                html.H6("MLST type: {}".format(mlst_type), className="table-header"),
            ])
        ]),
        html_test_tables(sample_data),
        res_div
    ])


def children_sample_list_report(dataframe):
    report = []
    row_index = 0
    for index, sample in \
            dataframe.iterrows():
        report.append(generate_sample_report(sample,
                                             row_index))
        row_index += 1
    return html.Div(report)

# callback
def samples_next_page(prev_ts, prev_ts2, next_ts, next_ts2, page_n, max_page):
    page_n = int(page_n)
    max_page = int(max_page)
    if max(prev_ts, prev_ts2) > max(next_ts, next_ts2):
        return str(max(page_n - 1, 0))
    elif max(next_ts, next_ts2) > max(prev_ts, prev_ts2):
        return str(min(page_n + 1, max_page))
    else:
        return '0'
