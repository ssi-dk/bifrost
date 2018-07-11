import dash_html_components as html
import dash_core_components as dcc
from components.images import list_of_images, get_species_color
from components.table import html_table, html_td_percentage
import plotly.graph_objs as go
import pandas as pd


def generate_sample_report(dataframe, sample, data_content, plot_data):
    return (
        html.Div(
            [
                html.A(id="sample-" + sample["name"]),
                html.H5(
                    sample["name"],
                    className="box-title"
                ),
                html_sample_tables(sample, data_content, className="row"),

                graph_sample_depth_plot(
                    sample, dataframe[dataframe.qcquickie_name_classified_species_1 == sample["qcquickie_name_classified_species_1"]]),
                graph_sample_cov_plot(sample, plot_data[data_content]["contig_coverage"])
            ],
            className="border-box"
        )
    )


def html_species_report(dataframe, species, data_content, plot_data, **kwargs):
    report = []
    for index, sample in dataframe.loc[dataframe["qcquickie_name_classified_species_1"] == species].iterrows():
        report.append(generate_sample_report(dataframe, sample, data_content, plot_data[sample["_id"]]))
    return html.Div(report, **kwargs)


def html_organisms_table(sample_data, **kwargs):
    percentages = [
        sample_data["qcquickie_percent_classified_species_1"],
        sample_data["qcquickie_percent_classified_species_2"],
        sample_data["qcquickie_percent_unclassified"]
    ]

    color_1 = get_species_color(
        sample_data["qcquickie_name_classified_species_1"])  # Default
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
                    html.I(sample_data["qcquickie_name_classified_species_1"])),
                html_td_percentage(percentages[0], color_1)
            ]),
            html.Tr([
                html.Td(
                    html.I(sample_data["qcquickie_name_classified_species_2"])),
                html_td_percentage(percentages[1], color_2)
            ]),
            html.Tr([
                html.Td("Unclassified"),
                html_td_percentage(percentages[2], color_u)
            ])
        ])
    ], **kwargs)


def html_sample_tables(sample_data, data_content, **kwargs):
    """Generate the tables for each sample containing submitter information,
       detected organisms etc. """
    genus = str(sample_data["qcquickie_name_classified_species_1"]).split()[
        0].lower()
    if "{}.svg".format(genus) in list_of_images:
        img = html.Img(src="/static/" + str(sample_data["qcquickie_name_classified_species_1"]).split()
                       [0].lower() + ".svg", className="svg_bact")
    else:
        img = []
    if type(sample_data["emails"]) is str:
        n_emails = len(sample_data["emails"].split(";"))
        if (n_emails > 1):
            emails = ", ".join(sample_data["emails"].split(";")[:2])
            if (n_emails > 2):
                emails += ", ..."
        else:
            emails = sample_data["emails"]
    else:
        emails = ''

    if data_content == "qcquickie":
        title = "QCQuickie Results"
        report = [
            html.Div([
                html_table([
                    [
                        "Number of contigs",
                        "{:,}".format(
                            sample_data["qcquickie_bin_contigs_at_1x"])
                    ],
                    [
                        "N50",
                        "{:,}".format(sample_data["qcquickie_N50"])
                    ],
                    [
                        "N75",
                        "{:,}".format(sample_data["qcquickie_N75"])
                    ],
                    [
                        "bin length at 1x depth",
                        "{:,}".format(
                            sample_data["qcquickie_bin_length_at_1x"])
                    ],
                    [
                        "bin length at 10x depth",
                        "{:,}".format(
                            sample_data["qcquickie_bin_length_at_10x"])
                    ],
                    [
                        "bin length at 25x depth",
                        "{:,}".format(
                            sample_data["qcquickie_bin_length_at_25x"])
                    ]
                ])
            ], className="six columns"),
            html_organisms_table(sample_data, className="six columns")
        ]
    elif data_content == "assembly":
        title = "Assembly Results"
        report = [
            html.Div([
                html_table([
                    [
                        "Number of contigs",
                        "{:,}".format(
                            sample_data["assembly_bin_contigs_at_1x"])
                    ],
                    [
                        "N50",
                        "{:,}".format(sample_data["assembly_N50"])
                    ],
                    [
                        "N75",
                        "{:,}".format(sample_data["assembly_N75"])
                    ],
                    [
                        "bin length at 1x depth",
                        "{:,}".format(
                            sample_data["assembly_bin_length_at_1x"])
                    ],
                    [
                        "bin length at 10x depth",
                        "{:,}".format(
                            sample_data["assembly_bin_length_at_10x"])
                    ],
                    [
                        "bin length at 25x depth",
                        "{:,}".format(
                            sample_data["assembly_bin_length_at_25x"])
                    ]
                ])
            ], className="six columns"),
            html_organisms_table(sample_data, className="six columns")
        ]
    else:
        title = "No report selected"
        report = []

    return html.Div([
        html.Div(img, className="bact_div"),
        html.H6("Run folder: " + sample_data["run_name"]),
        html.H6("Setup time: " + str(sample_data["setup_time"])),
        html.H5("Sample Sheet", className="table-header"),
        html.Div([
            html.Div([
                html_table([
                    ["Supplied name", sample_data["supplied_name"]],
                    ["Supplying lab", sample_data["supplying_lab"]],
                    ["Submitter emails", emails],
                    ["Provided species", html.I(
                        sample_data["provided_species"])]
                ])
            ], className="six columns"),
            html.Div([
                html.H6("User Comments", className="table-header"),
                sample_data["comments"]
            ], className="six columns"),
        ], className="row"),
        html.H5(title, className="table-header"),
        html.Div(report, className="row")
    ], **kwargs)


def graph_sample_depth_plot(sample, background_dataframe):
    # With real data, we should be getting sample data (where to put 1, 10
    # and 25x annotation) and the info for the rest of that species box.
    return dcc.Graph(
        id="coverage-1-" + sample['name'],
        figure={
            "data": [
                go.Box(
                    x=background_dataframe["qcquickie_bin_length_at_1x"],
                    text=background_dataframe["name"],
                    boxpoints="all",
                    jitter=0.3,
                    pointpos=-1.8,
                    marker=dict(color=get_species_color(
                        sample['qcquickie_name_classified_species_1']))
                )
                #{"x": [1, 2, 3], "y": [2, 4, 5], "type": "bar", "name": u"Montréal"},
            ],
            "layout": go.Layout(
                title="{}: Binned Depth 1x size".format(sample['name']),
                margin=go.Margin(
                    l=75,
                    r=50,
                    b=25,
                    t=50
                ),
                annotations=go.Annotations([
                    go.Annotation(
                        x=sample["qcquickie_bin_length_at_1x"],
                        y=0,
                        text="1x",
                        showarrow=True,
                        ax=40,
                        ay=0
                    ),
                    go.Annotation(
                        x=sample["qcquickie_bin_length_at_10x"],
                        y=0.0,
                        text="10x",
                        showarrow=True,
                        ax=0,
                        ay=40
                    ),
                    go.Annotation(
                        x=sample["qcquickie_bin_length_at_25x"],
                        y=0,
                        text="25x",
                        showarrow=True,
                        ax=0,
                        ay=-40
                    ),
                ])
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
                go.Scattergl(
                    x= df.total_length,
                    y= df.coverage,
                    text= df.index,
                    mode= 'markers'
                )
            ],
            "layout": go.Layout(
                title="{}: Binned Depth 1x size".format("something"),
                xaxis=dict(
                    type='log',
                    autorange=True
                ),
                yaxis=dict(
                    type='log',
                    autorange=True
                )
            )
        },
    )


def children_sample_list_report(filtered_df, data_content, plot_data):
    report = []
    for species in filtered_df.qcquickie_name_classified_species_1.unique():
        report.append(html.Div([
            html.A(id="species-cat-" + str(species).replace(" ", "-")),
            html.H4(html.I(str(species))),
            
            html_species_report(filtered_df, species, data_content, plot_data)
        ]))
    return report


