import dash_html_components as html
import dash_core_components as dcc
import dash_table as dt
from components.images import list_of_images
from components.table import html_table, html_td_percentage
import components.import_data as import_data
import components.global_vars as global_vars
import components.admin as admin
import dash_bootstrap_components as dbc
import pandas as pd
import math
import json

SAMPLE_PAGESIZE = 25

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
    test_path = "properties.stamper.summary." + test_name
    if test_path not in sample or pd.isnull(sample[test_path]):
        return "" # show nothing
        #return "test-missing"
    if sample.get(test_path, "").startswith("pass"):
        return "test-pass"
    elif sample.get(test_path, "").startswith("fail"):
        return "test-fail"
    else:
        return "test-warning"

def get_species_img(sample_data):
    genus = str(sample_data.get("properties.species_detection.summary.name_classified_species_1")).split()[
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
                        html.I(sample["properties.species_detection.summary.species"])
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
        sample_data.get("properties.species_detection.summary.percent_classified_species_1", math.nan),
        sample_data.get(
            "properties.species_detection.summary.percent_classified_species_2", math.nan),
        sample_data.get("properties.species_detection.summary.percent_unclassified", math.nan)
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
                    [html.I(sample_data.get("properties.species_detection.summary.name_classified_species_1", "No data")), " + Unclassified"], className="cell"),
                html_td_percentage(percentages[0] + percentages[2], color_0)
            ], className=check_test("whats_my_species:minspecies", sample_data) + " trow"),
            html.Tr([
                html.Td(
                    html.I(sample_data.get("properties.species_detection.summary.name_classified_species_1", "No data")), className="cell"),
                html_td_percentage(percentages[0], color_1)
            ], className="trow"),
            html.Tr([
                html.Td(
                    html.I(sample_data.get("properties.species_detection.summary.name_classified_species_2", "No data")), className="cell"),
                html_td_percentage(percentages[1], color_2)
            ], className="trow"),
            html.Tr([
                html.Td("Unclassified", className="cell"),
                html_td_percentage(percentages[2], color_u)
            ], className=check_test("whats_my_species:maxunclassified", sample_data) + " trow")
        ], className="bifrost-table")
    ], **kwargs)

def html_test_tables(sample_data):
    stamps_to_check = ["ssi_stamper", "supplying_lab_check"]
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

    stamp_rows = []
    for stamp in stamps_to_check:
        stamp_key = "stamps.{}.value".format(stamp)
        if stamp_key in sample_data and not pd.isnull(sample_data[stamp_key]):
            if str(sample_data[stamp_key]).startswith("pass"):
                stamp_class = "test-pass"
            elif str(sample_data[stamp_key]).startswith("fail"):
                stamp_class = "test-fail"
            else:
                stamp_class = ""
            stamp_rows.append({
                "list": [stamp, sample_data[stamp_key]],
                "className": stamp_class
            })

    if len(rows):
        test_table = html_table(rows)
    else:
        test_table = html.P("No failed tests.")
    
    if len(stamp_rows):
        stamp_table = html_table(stamp_rows)
    else:
        stamp_table = html.P("No stamps.")

    return dbc.Row([
        dbc.Col([
            html.H6("QC stamps", className="table-header"),
            stamp_table
        ]),
        dbc.Col([
            html.H6("Failed QC tests", className="table-header"),
            test_table
        ]),
    ])


def html_sample_tables(sample_data, **kwargs):
    """Generate the tables for each sample containing submitter information,
       detected organisms etc. """

    if "sample_sheet.sample_name" in sample_data:
        if "sample_sheet.emails" in sample_data and type(sample_data["sample_sheet.emails"]) is str:
            n_emails = len(sample_data["sample_sheet.emails"].split(";"))
            if (n_emails > 1):
                emails = ", ".join(
                    sample_data["sample_sheet.emails"].split(";")[:2])
                if (n_emails > 2):
                    emails += ", ..."
            else:
                emails = sample_data["sample_sheet.emails"]
        else:
            emails = ""
        sample_sheet_table = html_table([
                    ["Supplied name", sample_data.get("sample_sheet.sample_name","")],
                    ["User Comments", sample_data.get("sample_sheet.Comments","")],
                    ["Supplying lab", sample_data.get("sample_sheet.group", "")],
                    ["Submitter emails", emails],
                    {
                        "list": ["Provided species", html.I(
                            sample_data.get("sample_sheet.provided_species"))],
                        "className": check_test("properties.species_detection.summary:detectedspeciesmismatch", sample_data)
                    },
                    {
                        "list": ["Read file", 
                            str(sample_data["reads.R1"]).split("/")[-1]],
                        "className": check_test("base:readspresent", sample_data)
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
                        sample_data.get("properties.denovo_assembly.summary.filtered_reads_num", math.nan))
                ],
                "className": check_test("assemblatron:numreads", sample_data)
            },
            [
                "Number of contigs (1x cov.)",
                "{:,.0f}".format(
                    sample_data.get("properties.denovo_assembly.summary.bin_contigs_at_1x", math.nan))
            ],
            [
                "Number of contigs (10x cov.)",
                "{:,.0f}".format(
                    sample_data.get("properties.denovo_assembly.summary.bin_contigs_at_10x", math.nan))
            ],
            [
                "N50",
                "{:,}".format(sample_data.get("properties.denovo_assembly.summary.N50", math.nan))
            ],
            {
                "list": [
                    "Average coverage (1x)",
                    "{:,.2f}".format(
                        sample_data.get("properties.denovo_assembly.summary.bin_coverage_at_1x", math.nan))
                ],
                "className": check_test("assemblatron:avgcoverage", sample_data)
            },
            {
                "list": [
                    "Genome size at 1x depth",
                    "{:,.0f}".format(
                        sample_data.get("properties.denovo_assembly.summary.bin_length_at_1x", math.nan))
                ],
                "className": check_test("assemblatron:1xgenomesize", sample_data)
            },
            {
                "list": [
                    "Genome size at 10x depth",
                    "{:,.0f}".format(
                        sample_data.get("properties.denovo_assembly.summary.bin_length_at_10x", math.nan))
                ],
                "className": check_test("assemblatron:10xgenomesize", sample_data)
            },
            {
                "list": [
                    "Genome size 1x - 10x diff",
                    "{:,.0f}".format(
                        sample_data.get(
                            "properties.denovo_assembly.summary.bin_length_at_1x", math.nan)
                        - sample_data.get("properties.denovo_assembly.summary.bin_length_at_10x", math.nan)
                    )
                ],
                "className": check_test("assemblatron:1x10xsizediff", sample_data)
            },
            [
                "Genome size at 25x depth",
                "{:,.0f}".format(
                    sample_data.get("properties.denovo_assembly.summary.bin_length_at_25x", math.nan))
            ],
            [
                "Ambiguous sites",
                "{:,.0f}".format(
                    sample_data.get("properties.denovo_assembly.summary.snp_filter_10x_10%", math.nan))
            ] 
        ])
    ])
    expected_results = global_vars.expected_results
    any_results = False
    results = []
    for entry in expected_results:
        if "report.{}.data".format(entry) in sample_data:
            any_results = True
            results.append(html.Div([
                html.H6(sample_data["report.{}.title".format(entry)],
                        className="table-header"),
                html.Div(
                    dt.DataTable(
                        style_table={
                            'overflowX': 'scroll',
                            'overflowY': 'scroll',
                            'maxHeight': '480'
                        },

                        columns=sample_data["report.{}.columns".format(entry)],
                        data=sample_data["report.{}.data".format(entry)],
                        page_action='none'
                    ), className="grey-border")
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
    if "properties.mlst.summary.strain" in sample_data and not pd.isna(sample_data["properties.mlst.summary.strain"]):
        mlst_type = ", ".join(sample_data["properties.mlst.summary.strain"])


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
