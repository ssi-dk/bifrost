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
        html.H3("Sample Report"),
        html.Span("0", style={"display": "none"}, id="page-n"),
        
        html.Span(str(sample_n // SAMPLE_PAGESIZE),
                  style={"display": "none"}, id="max-page"),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Button(
                            "Previous page", id="prevpage", n_clicks_timestamp=0),
                    ]
                ),
                dbc.Col(
                    [
                        # html.H4("Page {} of {}".format(page_n + 1, max_page + 1))
                    ]
                ),
                dbc.Col(
                    [
                        html.Button(
                            "Next page", id="nextpage", n_clicks_timestamp=0),
                    ]
                ),
            ]
        ),

        html.Div(id="sample-report"),

        html.Div(dbc.Row(
            [
                dbc.Col(
                    [
                        html.Button(
                            "Previous page", id="prevpage2", n_clicks_timestamp=0),
                    ]
                ),
                dbc.Col(
                    [
                        # html.H4("Page {} of {}".format(page_n + 1, max_page + 1))
                    ]
                ),
                dbc.Col(
                    [
                        html.Button(
                            "Next page", id="nextpage2", n_clicks_timestamp=0),
                    ]
                ),
            ]
        ), className="mb-4"),
    ]

def check_test(test_name, sample):
    test_path = "sample_components.ssi_stamper.summary." + test_name
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
    genus = str(sample_data.get("sample_components.whats_my_species.summary.name_classified_species_1")).split()[
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
    img_div = dbc.CardTitle([
        html.Div(html.I(sample["properties.species"]),
                 className="d-inline-block mr-4"),
        img
    ], className="float-right")
    return (
        dbc.Card(
            [
                dbc.CardHeader([
                    html.A(id="sample-" + str(sample["name"])),
                    dbc.CardTitle(
                        sample["name"],
                        className="d-inline-block float-left"
                    ),
                    img_div,
                ]),
                dbc.CardBody([
                    html_sample_tables(sample),
                    admin.sample_radio_feedback(sample, n_sample)
                ])
                
            ]
        )
    )


def html_organisms_table(sample_data, **kwargs):
    percentages = [
        sample_data.get("sample_components.whats_my_species.summary.percent_classified_species_1", math.nan),
        sample_data.get(
            "sample_components.whats_my_species.summary.percent_classified_species_2", math.nan),
        sample_data.get("sample_components.whats_my_species.summary.percent_unclassified", math.nan)
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
                    [html.I(sample_data.get("sample_components.whats_my_species.summary.name_classified_species_1", "No data")), " + Unclassified"], className="cell"),
                html_td_percentage(percentages[0] + percentages[2], color_0)
            ], className=check_test("whats_my_species:minspecies", sample_data) + " trow"),
            html.Tr([
                html.Td(
                    html.I(sample_data.get("sample_components.whats_my_species.summary.name_classified_species_1", "No data")), className="cell"),
                html_td_percentage(percentages[0], color_1)
            ], className="trow"),
            html.Tr([
                html.Td(
                    html.I(sample_data.get("sample_components.whats_my_species.summary.name_classified_species_2", "No data")), className="cell"),
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
        if key.startswith("sample_components.ssi_stamper.summary.whats_my_species") \
                or key.startswith("sample_components.ssi_stamper.summary.assemblatron"):
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
                        "className": check_test("sample_components.whats_my_species.summary:detectedspeciesmismatch", sample_data)
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
                        sample_data.get("sample_components.assemblatron.summary.filtered_reads_num", math.nan))
                ],
                "className": check_test("assemblatron:numreads", sample_data)
            },
            [
                "Number of contigs (1x cov.)",
                "{:,.0f}".format(
                    sample_data.get("sample_components.assemblatron.summary.bin_contigs_at_1x", math.nan))
            ],
            [
                "Number of contigs (10x cov.)",
                "{:,.0f}".format(
                    sample_data.get("sample_components.assemblatron.summary.bin_contigs_at_10x", math.nan))
            ],
            [
                "N50",
                "{:,}".format(sample_data.get("sample_components.assemblatron.summary.N50", math.nan))
            ],
            {
                "list": [
                    "Average coverage (1x)",
                    "{:,.2f}".format(
                        sample_data.get("sample_components.assemblatron.summary.bin_coverage_at_1x", math.nan))
                ],
                "className": check_test("assemblatron:avgcoverage", sample_data)
            },
            {
                "list": [
                    "Genome size at 1x depth",
                    "{:,.0f}".format(
                        sample_data.get("sample_components.assemblatron.summary.bin_length_at_1x", math.nan))
                ],
                "className": check_test("assemblatron:1xgenomesize", sample_data)
            },
            {
                "list": [
                    "Genome size at 10x depth",
                    "{:,.0f}".format(
                        sample_data.get("sample_components.assemblatron.summary.bin_length_at_10x", math.nan))
                ],
                "className": check_test("assemblatron:10xgenomesize", sample_data)
            },
            {
                "list": [
                    "Genome size 1x - 10x diff",
                    "{:,.0f}".format(
                        sample_data.get(
                            "sample_components.assemblatron.summary.bin_length_at_1x", math.nan)
                        - sample_data.get("sample_components.assemblatron.summary.bin_length_at_10x", math.nan)
                    )
                ],
                "className": check_test("assemblatron:1x10xsizediff", sample_data)
            },
            [
                "Genome size at 25x depth",
                "{:,.0f}".format(
                    sample_data.get("sample_components.assemblatron.summary.bin_length_at_25x", math.nan))
            ],
            [
                "Ambiguous sites",
                "{:,.0f}".format(
                    sample_data.get("sample_components.assemblatron.summary.snp_filter_10x_10%", math.nan))
            ] 
        ])
    ])
    resresults = False
    resfinder = sample_data.get('sample_components.ariba_resfinder.summary.ariba_resfinder', [])
    if type(resfinder) == list and len(resfinder):
        resresults = True
        resfinder_div = html.Div(
            dt.DataTable(
                style_table={
                    'overflowX': 'scroll',
                    'overflowY': 'scroll',
                    'maxHeight': '480'
                },

                columns=global_vars.finder_columns,
                data=resfinder,
                pagination_mode=False
            ), className="grey-border")
    elif (sample_data.get("sample_components.ariba_resfinder.status", "") == "Success" and
         (type(resfinder) == float or resfinder is None or not len(resfinder))):
        resfinder_div = html.P("No antibiotic resistance genes found")
    else:
        resfinder_div = html.P("Resfinder not run")
    
    plasmidfinder = sample_data.get('ariba_plasmidfinder.ariba_plasmidfinder', [])
    if type(plasmidfinder) == list and len(plasmidfinder):
        resresults = True
        plasmidfinder_div = html.Div(
            dt.DataTable(
                style_table={
                    'overflowX': 'scroll',
                    'overflowY': 'scroll',
                    'maxHeight': '480'
                },
                columns=global_vars.finder_columns,
                data=plasmidfinder,
                pagination_mode=False
            ), className="grey-border")
    elif (sample_data.get("sample_components.ariba_plasmidfinder.status", "") == "Success" and
         (type(plasmidfinder) == float or
          plasmidfinder is None or
          not len(plasmidfinder))):
        plasmidfinder_div = html.P("No replicons found")
    else:
        plasmidfinder_div = html.P("Plasmidfinder not run")

    virulencefinder = sample_data.get('ariba_virulencefinder.ariba_virulencefinder', [])
    if type(virulencefinder) == list and len(virulencefinder):
        resresults = True
        virulencefinder_div = html.Div(
            dt.DataTable(
                style_table={
                    'overflowX': 'scroll',
                    'overflowY': 'scroll',
                    'maxHeight': '480'
                },
                columns=global_vars.finder_columns,
                data=virulencefinder,
                pagination_mode=False
            ), className="grey-border")
    elif sample_data.get("sample_components.ariba_virulencefinder.status", "") == "Success" and (type(virulencefinder) == float or virulencefinder is None or not len(virulencefinder)):
        virulencefinder_div = html.P("No virulence markers found")
    else:
        virulencefinder_div = html.P("Virulencefinder not run")

    mlst_data = sample_data.get(
        'sample_components.ariba_mlst.summary.mlst_report', "")
    if type(mlst_data) == str and len(mlst_data):
        resresults = True
        mlst_dict = {}
        mlst_fields = mlst_data.split(",")
        for field in mlst_fields:
            key, value = field.split(":")
            mlst_dict[key] = value

        columns = [{"name": i, "id": i} for i in mlst_dict.keys()]
        mlst_div = html.Div(
            dt.DataTable(
                style_table={
                    'overflowX': 'scroll',
                    'overflowY': 'scroll',
                    'maxHeight': '480'
                },
                columns=columns,
                data=[mlst_dict],
                pagination_mode=False
            ), className="grey-border")
    else:
        mlst_div = html.P("MLST not run")

    mlst_db = sample_data.get("sample_components.ariba_mlst.summary.mlst_db", "")

    # Replace with the ariba_res, ariba_plas and ariba_vir when migrating to them
    if (sample_data.get("sample_components.ariba_resfinder.status", "") == "Success" or
        sample_data.get("sample_components.ariba_plasmidfinder.status", "") == "Success" or
        sample_data.get("sample_components.ariba_mlst.status", "") == "Success" or
        sample_data.get("sample_components.ariba_virulencefinder.status", "") == "Success"):
        res_analysis_not_run = False
    else:
        res_analysis_not_run = True

    if resresults:
        res_div = html.Details([
            html.Summary("ResFinder/PlasmidFinder/VirulenceFinder/MLST (click to show)"),
            html.Div([
                html.Div([
                    html.Div([
                        html.H6("ResFinder", className="table-header"),
                        resfinder_div
                    ], className="six columns"),
                    html.Div([
                        html.H6("VirulenceFinder", className="table-header"),
                        virulencefinder_div
                    ], className="six columns")
                ], className="row"),
                html.Div([
                    html.Div([
                        html.H6("PlasmidFinder", className="table-header"),
                        plasmidfinder_div
                    ], className="six columns"),
                    html.Div([
                        html.H6("MLST ({})".format(mlst_db), className="table-header"),
                        mlst_div
                    ], className="six columns")
                ], className="row")
            ])
        ])
    elif not resresults and not res_analysis_not_run:
        res_div = html.Div(html.P("No antibiotic resistances, replicons or virulence markers found for this sample."))
    else:
        res_div = html.Div(
            html.P("Resfinder, plasmidfinder and virulencefinder were not run."))

    mlst_type = "ND"
    if "sample_components.ariba_mlst.summary.mlst_report" in sample_data and sample_data["sample_components.ariba_mlst.summary.mlst_report"] is not None:
        mlst_report_string = sample_data["sample_components.ariba_mlst.summary.mlst_report"]
        if "," in str(mlst_report_string):
            mlst_text_split = mlst_report_string.split(",", 1)
            mlst_type = mlst_text_split[0].split(":",1)[1]


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

def generate_sample_folder(sample_ids):
    """Generates a script string """
    samples = import_data.get_samples(sample_ids)
    # Access samples by their id
    samples_by_ids = { str(s["_id"]) : s for s in samples }
    assemblies = import_data.get_assemblies_paths(sample_ids)
    reads_script = "mkdir samples\ncd samples\n"
    assemblies_script = "mkdir assemblies\ncd assemblies\n"
    reads_errors = []
    assemblies_errors = []
    for assembly in assemblies:
        sample = samples_by_ids[str(assembly["sample"]["_id"])]
        try:
            assemblies_script += "#{}\nln -s {} {}\n".format(
                sample["name"],
                assembly["path"] + "/contigs.fasta",
                sample["name"] + "_contigs.fasta")
        except KeyError as e:
            assemblies_errors.append("Missing data for sample: {} - {}. In database:\n{}".format(
                sample.get("name", "None"),
                assembly["_id"],
                assembly.get("path", "No data")
            ))
    
    if len(assemblies_errors):
        assemblies_html = [
            html.H5(
                "Use this script to generate a folder with all the assemblies linked in it."),
            "A few errors occurred locating the contigs. If you need more info, " +
            "please contact an admin.",
            html.Pre("\n".join(assemblies_errors), className="error-pre"),
            html.Pre(assemblies_script, className="folder-pre")
        ]
    else:
        assemblies_html = [
            html.H5(
                "Use this script to generate a folder with all the assemblies linked in it."),
            html.Pre(assemblies_script, className="folder-pre")
        ]

    for sample in samples:
        try:
            reads_script += "#{}\nln -s {} .\nln -s {} .\n".format(
                sample["name"],
                sample["reads"]["R1"],
                sample["reads"]["R2"])
        except KeyError as e:
            reads_errors.append("Missing data for sample: {} - {}. In database:\n{}".format(
                sample.get("name", "None"),
                sample["_id"],
                sample.get("reads", "No data")
                ))
    if len(reads_errors):
        reads_html = [
            html.H5("Use this script to generate a folder with all the sample reads linked in it."),
            "A few errors occurred locating the read paths. If you need more info, " +
            "please contact an admin.",
            html.Pre("\n".join(reads_errors), className="error-pre"),
            html.Pre(reads_script, className="folder-pre")
        ]
    else:
        reads_html = [
            html.H5("Use this script to generate a folder with all the sample reads linked in it."),
            html.Pre(reads_script, className="folder-pre")
            ]
    return html.Div([html.Div(assemblies_html), html.Div(reads_html)])

