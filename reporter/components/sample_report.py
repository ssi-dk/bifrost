import dash_html_components as html
import dash_core_components as dcc
import dash_table as dt
from components.images import list_of_images, get_species_color
from components.table import html_table, html_td_percentage
from import_data import get_read_paths
import plotly.graph_objs as go
import pandas as pd
import math
import json

def check_test(test_name, sample):
    test_path = "ssi_stamper." + test_name
    if test_path not in sample or pd.isnull(sample[test_path]):
        return "" # show nothing
        #return "test-missing"
    if sample[test_path].startswith("pass"):
        return "test-pass"
    elif sample[test_path].startswith("fail"):
        return "test-fail"
    else:
        return "test-warning"

def get_species_img(sample_data):
    genus = str(sample_data.get("whats_my_species.name_classified_species_1")).split()[
        0].lower()
    if "{}.svg".format(genus) in list_of_images:
        img = html.Img(
            src="/assets/img/" + genus + ".svg",
            className="svg_bact"
        )
    else:
        img = None
    return img

def generate_sample_report(dataframe, sample, background):
    plot = graph_sample_depth_plot(
        sample,
        dataframe[dataframe["species"]
                == sample["species"]],
        background
    )
    tests = html_test_table(sample, className="row")

    img = get_species_img(sample)
    if img is not None:
        img_div = html.Div(img, className="box-title bact grey-border")
    else:
        img_div = None
    return (
        html.Div(
            [
                html.A(id="sample-" + sample["name"]),
                html.H5(
                    sample["name"],
                    className="box-title"
                ),
                img_div,
                html_sample_tables(sample, className="row"),
                # plot,
                html.H6("Failed QC Tests", className="table-header"),
                tests
            ],
            className="border-box"
        )
    )


def html_species_report(dataframe, species, species_plot_data, **kwargs):
    report = []
    for index, sample in \
      dataframe.loc[dataframe["species"] == species].iterrows():
        report.append(generate_sample_report(dataframe,
                                             sample,
                                             species_plot_data))
    return html.Div(report, **kwargs)


def html_organisms_table(sample_data, **kwargs):
    percentages = [
        sample_data.get("whats_my_species.percent_classified_species_1", math.nan),
        sample_data.get(
            "whats_my_species.percent_classified_species_2", math.nan),
        sample_data.get("whats_my_species.percent_unclassified", math.nan)
    ]

    color_0 = "#b3ccc1"

    color_1 = get_species_color(
        sample_data.get("whats_my_species.name_classified_species_1"))  # Default
#    color_2 = COLOR_DICT.get(
#        sample_data["name_classified_species_2"], "#f3bbd3")  # Default
    color_2 = "#f3bbd3"  # Default

#   color_u = COLOR_DICT.get("", "#fee1cd")  # Default
    color_u = "#fee1cd"  # Default

    return html.Div([
        html.Table([
            html.Tr([
                html.Td(
                    [html.I(sample_data.get("whats_my_species.name_classified_species_1", "No data")), " + Unclassified"], className="cell"),
                html_td_percentage(percentages[0] + percentages[2], color_0)
            ], className=check_test("whats_my_species:minspecies", sample_data) + " trow"),
            html.Tr([
                html.Td(
                    html.I(sample_data.get("whats_my_species.name_classified_species_1", "No data")), className="cell"),
                html_td_percentage(percentages[0], color_1)
            ], className="trow"),
            html.Tr([
                html.Td(
                    html.I(sample_data.get("whats_my_species.name_classified_species_2", "No data")), className="cell"),
                html_td_percentage(percentages[1], color_2)
            ], className="trow"),
            html.Tr([
                html.Td("Unclassified", className="cell"),
                html_td_percentage(percentages[2], color_u)
            ], className=check_test("whats_my_species:maxunclassified", sample_data) + " trow")
        ])
    ], **kwargs)

def html_test_table(sample_data,**kwargs):
    rows = []
    for key, value in sample_data.items():
        if key.startswith("ssi_stamper.whats_my_species") \
        or key.startswith("ssi_stamper.assemblatron"):
            if pd.isnull(value):
                value = "nan"
            name_v = key.split(":")
            values = value.split(':')
            if len(values) == 3:
                if values[0] != "pass":
                    rows.append([name_v[1].capitalize(), values[1]])
            if (key.endswith(".action")):
                rows.append(["QC Action", value])
    if len(rows):
        return html.Div(html_table(rows, className="twelve columns"), **kwargs)
    else:
        return html.Div(html.P("No failed tests."), **kwargs)

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
                    ["Supplied name", sample_data["sample_sheet.sample_name"]],
                    ["User Comments", sample_data["sample_sheet.Comments"]],
                    ["Supplying lab", sample_data["sample_sheet.group"]],
                    ["Submitter emails", emails],
                    {
                        "list": ["Provided species", html.I(
                            sample_data["sample_sheet.provided_species"])],
                        "className": check_test("whats_my_species:detectedspeciesmismatch", sample_data)
                    },
                    {
                        "list": ["Read file", 
                            str(sample_data["R1"]).split("/")[-1]],
                        "className": check_test("base:readspresent", sample_data)
                    }
                ])
    else:
        sample_sheet_table = []

    title = "Assemblatron Results"
    table = html.Div([
        html_table([
            {
                "list": [
                    "Number of filtered reads",
                    "{:,.0f}".format(
                        sample_data.get("assemblatron.filtered_reads_num", math.nan))
                ],
                "className": check_test("assemblatron:numreads", sample_data)
            },
            [
                "Number of contigs (1x cov.)",
                "{:,.0f}".format(
                    sample_data.get("assemblatron.bin_contigs_at_1x", math.nan))
            ],
            [
                "Number of contigs (10x cov.)",
                "{:,.0f}".format(
                    sample_data.get("assemblatron.bin_contigs_at_10x", math.nan))
            ],
            [
                "N50",
                "{:,}".format(sample_data.get("assemblatron.N50", math.nan))
            ],
            {
                "list": [
                    "Average coverage (1x)",
                    "{:,.2f}".format(
                        sample_data.get("assemblatron.bin_coverage_at_1x", math.nan))
                ],
                "className": check_test("assemblatron:avgcoverage", sample_data)
            },
            {
                "list": [
                    "Genome size at 1x depth",
                    "{:,.0f}".format(
                        sample_data.get("assemblatron.bin_length_at_1x", math.nan))
                ],
                "className": check_test("assemblatron:1xgenomesize", sample_data)
            },
            {
                "list": [
                    "Genome size at 10x depth",
                    "{:,.0f}".format(
                        sample_data.get("assemblatron.bin_length_at_10x", math.nan))
                ],
                "className": check_test("assemblatron:10xgenomesize", sample_data)
            },
            {
                "list": [
                    "Genome size 1x - 10x diff",
                    "{:,.0f}".format(
                        sample_data.get(
                            "assemblatron.bin_length_at_1x", math.nan)
                        - sample_data.get("assemblatron.bin_length_at_10x", math.nan)
                    )
                ],
                "className": check_test("assemblatron:1x10xsizediff", sample_data)
            },
            [
                "Genome size at 25x depth",
                "{:,.0f}".format(
                    sample_data.get("assemblatron.bin_length_at_25x", math.nan))
            ],
            [
                "Ambiguous sites",
                "{:,.0f}".format(
                    sample_data.get("assemblatron.snp_filter_10x_10%", math.nan))
            ] 
        ])
    ])

    resfinder = sample_data.get('analyzer.ariba_resfinder', "[{}]")
    if type(resfinder) == str:
        resfinder = resfinder.replace("'", "\"")
        resfinder = json.loads(resfinder)
        columns = []
        if len(resfinder):
            columns = [{"name": i, "id": i} for i in resfinder[0].keys()]
            resfinder_div = html.Div(
                dt.DataTable(
                    style_table={
                        'overflowX': 'scroll',
                        'overflowY': 'scroll',
                        'maxHeight': '480'
                    },

                    columns=columns,
                    data=resfinder,
                    pagination_mode=False
                ), className="grey-border")
    if type(resfinder) == float or not len(resfinder):
        resfinder_div = html.P("No antibiotic resistance genes found")
    
    return html.Div([
        html.Div([
            html.Div([
                html.H6("Sample Sheet", className="table-header"),
                sample_sheet_table,
                html.H6("Detected Organisms", className="table-header"),
                html_organisms_table(sample_data)
            ], className="six columns"),
            html.Div([
                html.H6(title, className="table-header"),
                table,
                html.H6("MLST", className="table-header"),
                html_table([[sample_data.get("analyzer.mlst_report", "No results")]])
            ], className="six columns")
        ], className="row"),
        html.Div([
            html.Div([
                html.H6("Resfinder", className="table-header"),
                resfinder_div
            ], className="twelve columns")
        ], className="row")
    ], **kwargs)


def graph_sample_depth_plot(sample, run_species, background):
    # With real data, we should be getting sample data (where to put 1, 10
    # and 25x annotation) and the info for the rest of that species box.
    return dcc.Graph(
        id="coverage-1-" + sample["_id"],
        figure={
            "data": [
                go.Box(
                    x=run_species.get("assemblatron.bin_length_at_1x"),
                    text=run_species["name"],
                    name="Current run",
                    showlegend=False,
                    boxpoints="all",
                    pointpos=-1.8,
                    jitter=0.3,
                    marker=dict(
                        size=4,
                        color=get_species_color(
                            sample.get("whats_my_species.name_classified_species_1"))
                    )
                ),
                go.Box(
                    x=background,
                    boxpoints="all",
                    showlegend=False,
                    name="Prev. runs",
                    jitter=0.3,
                    pointpos=-1.8,
                    marker=dict(
                        color="black",
                        size=4
                    )
                )
            ],
            "layout": go.Layout(
                title="{}: Genome size 1x depth".format(sample["name"]),
                hovermode="closest",
                margin=go.layout.Margin(
                    l=75,
                    r=50,
                    b=25,
                    t=50
                ),
                annotations=[
                    dict(
                        x=sample.get("assemblatron.bin_length_at_1x", []),
                        y=0,
                        text="1x",
                        showarrow=True,
                        ax=35,
                        ay=0
                    ),
                    dict(
                        x=sample.get("assemblatron.bin_length_at_10x", []),
                        y=0.0,
                        text="10x",
                        showarrow=True,
                        ax=0,
                        ay=35
                    ),
                    dict(
                        x=sample.get("assemblatron.bin_length_at_25x", []),
                        y=0,
                        text="25x",
                        showarrow=True,
                        ax=-35,
                        ay=0
                    ),
                ]
            )
        },
        style={"height": "200px"}
    )


def children_sample_list_report(filtered_df, plot_data):
    report = []
    for species in filtered_df["species"].unique():
        report.append(html.Div([
            html.A(id="species-cat-" + str(species).replace(" ", "-")),
            html.H4(html.I(str(species))),
            html_species_report(filtered_df, species,[])
            # html_species_report(filtered_df, species,
                                # plot_data.get(species, []))
        ]))
    return report

def generate_sample_folder(samples):
    """Generates a script string """
    reads = get_read_paths(samples)
    script = "mkdir samples\ncd samples\n"
    errors = []
    for sample in reads:
        try:
            script += "#{}\nln -s {} .\nln -s {} .\n".format(
                sample["name"],
                sample["reads"]["R1"],
                sample["reads"]["R2"])
        except KeyError as e:
            errors.append("Missing data for sample: {} - {}. In database:\n{}".format(
                sample.get("name", "None"),
                sample["_id"],
                sample.get("reads", "No data")
                ))
    if len(errors):
        return [
            html.H5("Use this script to generate a folder with all the sample reads linked in it."),
            "A few errors occurred locating the read paths. If you need more info, " +
            "please contact an admin.",
            html.Pre("\n".join(errors), className="error-pre"),
            html.Pre(script, className="folder-pre")
        ]
    else:
            
        return [
            html.H5("Use this script to generate a folder with all the sample reads linked in it."),
            html.Pre(script, className="folder-pre")
            ]


