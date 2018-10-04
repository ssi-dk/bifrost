import dash_html_components as html
import dash_table_experiments as dt
import dash_core_components as dcc
from components.images import list_of_images, get_species_color
from components.table import html_table, html_td_percentage
from import_data import get_read_paths
import plotly.graph_objs as go
import pandas as pd
import math

def check_test(test_name, sample):
    test_path = "testomatic." + test_name
    if test_path not in sample or pd.isnull(sample[test_path]):
        return "" # show nothing
        #return "test-missing"
    if sample[test_path].startswith("pass"):
        return "test-pass"
    elif sample[test_path].startswith("fail"):
        return "test-fail"
    else:
        return "test-warning"

def generate_sample_report(dataframe, sample, data_content, background):
    if data_content in ["qcquickie", "assemblatron"]:
        plot = graph_sample_depth_plot(
            sample,
            dataframe[dataframe["species"]
                    == sample["species"]],
            background
        )
        tests = html_test_table(sample, data_content, className="row")
    else:
        plot = None
        tests = None
    return (
        html.Div(
            [
                html.A(id="sample-" + sample["name"]),
                html.H5(
                    sample["name"],
                    className="box-title"
                ),
                html_sample_tables(sample, data_content, className="row"),

                plot,
                html.H6("QC Tests", className="table-header"),
                tests
            ],
            className="border-box"
        )
    )


def html_species_report(dataframe, species, data_content, species_plot_data, **kwargs):
    report = []
    for index, sample in \
      dataframe.loc[dataframe["species"] == species].iterrows():
        report.append(generate_sample_report(dataframe,
                                             sample,
                                             data_content,
                                             species_plot_data))
    return html.Div(report, **kwargs)


def html_organisms_table(sample_data, **kwargs):
    percentages = [
        sample_data.get("whats_my_species.percent_classified_species_1", math.nan),
        sample_data.get(
            "whats_my_species.percent_classified_species_2", math.nan),
        sample_data.get("whats_my_species.percent_unclassified", math.nan)
    ]

    color_1 = get_species_color(
        sample_data.get("whats_my_species.name_classified_species_1"))  # Default
#    color_2 = COLOR_DICT.get(
#        sample_data["name_classified_species_2"], "#f3bbd3")  # Default
    color_2 = "#f3bbd3"  # Default

#   color_u = COLOR_DICT.get("", "#fee1cd")  # Default
    color_u = "#fee1cd"  # Default

    return html.Div([
        html.H6("Detected Organisms", className="table-header"),
        html.Table([
            html.Tr([
                html.Td(
                    html.I(sample_data.get("whats_my_species.name_classified_species_1", "No data")), className="cell"),
                html_td_percentage(percentages[0], color_1)
            ], className=check_test("whats_my_species.minspecies", sample_data) + " trow"),
            html.Tr([
                html.Td(
                    html.I(sample_data.get("whats_my_species.name_classified_species_2", "No data")), className="cell"),
                html_td_percentage(percentages[1], color_2)
            ], className="trow"),
            html.Tr([
                html.Td("Unclassified", className="cell"),
                html_td_percentage(percentages[2], color_u)
            ], className=check_test("whats_my_species.maxunclassified", sample_data) + " trow")
        ])
    ], **kwargs)

def html_test_table(sample_data, data_content,**kwargs):
    rows = []
    for key, value in sample_data.items():
        if key.startswith("testomatic.whats_my_species") \
        or key.startswith("testomatic." + data_content):
            if pd.isnull(value):
                value = "nan"
            name_v = key.split(".")
            values = value.split(':')
            if len(values) == 3:
                if values[0] != "pass":
                    rows.append([name_v[2].capitalize(), "{}. Reason: {}. Detected value: {}".format(
                        values[0], values[1], values[2])])
            if (key.endswith(".action")):
                rows.append(["QC Action", value])
    return html.Div(html_table(rows, className="twelve columns"), **kwargs)

def html_sample_tables(sample_data, data_content, **kwargs):
    """Generate the tables for each sample containing submitter information,
       detected organisms etc. """
    genus = str(sample_data.get("whats_my_species.name_classified_species_1")).split()[
        0].lower()
    if "{}.svg".format(genus) in list_of_images:
        img = html.Img(
            src="/assets/img/" + genus + ".svg",
            className="svg_bact"
        )
    else:
        img = []

    if pd.isna(sample_data["runs"]):
        runs = "Run: None"
    else:
        runs = "Run: " + ",".join(sample_data["runs"])
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
        sample_sheet_div = [
            html.H5("Sample Sheet", className="table-header"),
            html.Div([
                html.Div([
                    html_table([
                        ["Supplied name", sample_data["sample_sheet.sample_name"]],
                        ["Supplying lab", sample_data["sample_sheet.group"]],
                        ["Submitter emails", emails],
                        {
                            "list": ["Provided species", html.I(
                                sample_data["sample_sheet.provided_species"])],
                            "className": check_test("whats_my_species.submitted==detected", sample_data)
                        }
                    ])
                ], className="six columns"),
                html.Div([
                    html.H6("User Comments", className="table-header"),
                    sample_data["sample_sheet.Comments"]
                ], className="six columns"),
            ], className="row"),
        ]
    else:
        sample_sheet_div = []

    if data_content == "qcquickie":
        title = "QCQuickie Results"
        report = [
            html.Div([
                html_table([
                    [
                        "Number of contigs",
                        "{:,}".format(
                            sample_data.get("qcquickie.bin_contigs_at_1x", math.nan))
                    ],
                    [
                        "N50",
                        "{:,}".format(
                            sample_data.get("qcquickie.N50", math.nan))
                    ],
                    {
                        "list": [
                            "bin length at 1x depth",
                            "{:,}".format(
                                sample_data.get("qcquickie.bin_length_at_1x", math.nan))
                        ],
                        "className": check_test("qcquickie.1xgenomesize", sample_data)
                    },
                    [
                        "average coverage (1x)",
                        "{:,.2f}".format(
                            sample_data.get("qcquickie.bin_coverage_at_1x", math.nan))
                    ],
                    [
                        "bin length at 25x depth",
                        "{:,}".format(
                            sample_data.get("qcquickie.bin_length_at_25x", math.nan))
                    ],
                    {
                        "list": [
                            "bin length 1x - 25x diff",
                            "{:,}".format(
                                sample_data.get("qcquickie.bin_length_at_1x", math.nan) \
                                - sample_data.get("qcquickie.bin_length_at_25x", math.nan)
                                )
                        ],
                        "className": check_test("qcquickie.1x25xsizediff", sample_data)
                    }
                ])
            ], className="six columns"),
            html_organisms_table(sample_data, className="six columns")
        ]
    elif data_content == "assemblatron":
        title = "Assemblatron Results"
        report = [
            html.Div([
                html_table([
                    [
                        "Number of contigs",
                        "{:,}".format(
                            sample_data.get("assemblatron.bin_contigs_at_1x", math.nan))
                    ],
                    [
                        "N50",
                        "{:,}".format(sample_data.get("assemblatron.N50", math.nan))
                    ],
                    {
                        "list": [
                            "bin length at 1x depth",
                            "{:,}".format(
                                sample_data.get("assemblatron.bin_length_at_1x", math.nan))
                        ],
                        "className": check_test("assemblatron.1xgenomesize", sample_data)
                    },
                    [
                        "average coverage (1x)",
                        "{:,.2f}".format(
                            sample_data.get("assemblatron.bin_coverage_at_1x", math.nan))
                    ],
                    [
                        "bin length at 10x depth",
                        "{:,}".format(
                            sample_data.get("assemblatron.bin_length_at_10x", math.nan))
                    ],
                    [
                        "bin length at 25x depth",
                        "{:,}".format(
                            sample_data.get("assemblatron.bin_length_at_25x", math.nan))
                    ],
                    {
                        "list": [
                            "bin length 1x - 25x diff",
                            "{:,}".format(
                                sample_data.get(
                                    "assemblatron.bin_length_at_1x", math.nan)
                                - sample_data.get("assemblatron.bin_length_at_25x", math.nan)
                            )
                        ],
                        "className": check_test("assemblatron.1x25xsizediff", sample_data)
                    }
                ])
            ], className="six columns"),
            html_organisms_table(sample_data, className="six columns")
        ]
    elif data_content == "analyzer":
        title = "Resfinder Results"
        resfinder = sample_data.get('analyzer.ariba_resfinder', None)
        if (isinstance(resfinder, list) and len(resfinder)):
            header = list(resfinder[0].keys())
            rows = [list(row.values()) for row in resfinder]
            report = [
                html.Div([
                    dt.DataTable(
                        rows=resfinder,
                        editable=False,
                        sortable=False,
                        column_widths=[100]*len(header)
                        ),
                    #html_table([header] + rows)
                ], className="twelve columns")
            ]
        else:
            report = ["No results found"]
    else:
        title = "No report selected"
        report = []

    return html.Div([
        html.Div(img, className="bact_div"),
        html.H6(runs),
        html.Div(sample_sheet_div),
        html.H5(title, className="table-header"),
        html.Div(report, className="row")
    ], **kwargs)


def graph_sample_depth_plot(sample, run_species, background):
    # With real data, we should be getting sample data (where to put 1, 10
    # and 25x annotation) and the info for the rest of that species box.
    return dcc.Graph(
        id="coverage-1-" + sample["_id"],
        figure={
            "data": [
                go.Box(
                    x=run_species.get("qcquickie.bin_length_at_1x"),
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
                title="{}: Binned Depth 1x size".format(sample["name"]),
                hovermode="closest",
                margin=go.layout.Margin(
                    l=75,
                    r=50,
                    b=25,
                    t=50
                ),
                annotations=[
                    dict(
                        x=sample.get("qcquickie.bin_length_at_1x", []),
                        y=0,
                        text="1x",
                        showarrow=True,
                        ax=35,
                        ay=0
                    ),
                    dict(
                        x=sample.get("qcquickie.bin_length_at_10x", []),
                        y=0.0,
                        text="10x",
                        showarrow=True,
                        ax=0,
                        ay=35
                    ),
                    dict(
                        x=sample.get("qcquickie.bin_length_at_25x", []),
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

def graph_sample_cov_plot(sample, sample_coverage):
    df = pd.DataFrame.from_dict(sample_coverage, orient="index")
    return dcc.Graph(
        id="something" + sample["name"],
        figure={
            "data": [
                go.Scatter(
                    x= df.total_length,
                    y= df.coverage,
                    text= df.index,
                    mode= "markers"
                )
            ],
            "layout": go.Layout(
                title="{}: Binned Depth 1x size".format("something"),
                xaxis=dict(
                    type="log",
                    autorange=True
                ),
                yaxis=dict(
                    type="log",
                    autorange=True
                )
            )
        },
    )


def children_sample_list_report(filtered_df, data_content, plot_data):
    report = []
    for species in filtered_df["species"].unique():
        report.append(html.Div([
            html.A(id="species-cat-" + str(species).replace(" ", "-")),
            html.H4(html.I(str(species))),
            html_species_report(filtered_df, species, data_content, plot_data.get(species,[]))
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
            "A few errors occurred locating the read paths. If you need more info, " +
            "please contact an admin.",
            html.Pre("\n".join(errors), style={
                "border": "1px solid red", "padding": "1em", "marginBottom": "20px"}),
            html.Pre(script, style={
                  "border": "1px solid black", "padding": "1em"})
        ]
    else:
        return [html.Pre(script, style={"border": "1px solid black", "padding": "1em"})]


