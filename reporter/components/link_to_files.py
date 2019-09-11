
import dash_html_components as html

import components.import_data as import_data


def link_to_files(sample_store):
    """Generates a script string """
    sample_ids = [str(s["_id"]) for s in sample_store]
    samples = import_data.filter_all(sample_ids=sample_ids,
                                      projection={"name": 1,
                                                  "properties.datafiles.summary.paired_reads": 1})
    samples["_id"] = samples["_id"].astype(str)
    samples = samples.set_index("_id")
    # Access samples by their id
    assemblies = import_data.get_assemblies_paths(sample_ids)
    reads_script = "mkdir samples\ncd samples\n"
    assemblies_script = "mkdir assemblies\ncd assemblies\n"
    reads_errors = []
    assemblies_errors = []
    for assembly in assemblies:
        sample = samples.loc[str(assembly["sample"]["_id"]), : ]
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
            html.Div([
                html.Div([
                    html.H6("Assemblies",
                            className="m-0 font-weight-bold text-primary"),
                ], className="card-header py-3"),
                html.Div([
                    "A few errors occurred locating the contigs. If you need more info, " +
                    "please contact an admin.",
                    html.Pre("\n".join(assemblies_errors),
                             className="error-pre"),
                    html.Pre(assemblies_script, className="folder-pre")
                ], className="card-body")
            ], className="card mb-4 shadow")
        ]
    else:
        assemblies_html = [
            html.Div([
                html.Div([
                    html.H6("Assemblies",
                            className="m-0 font-weight-bold text-primary"),
                ], className="card-header py-3"),
                html.Div([
                    html.Pre(assemblies_script, className="folder-pre")
                ], className="card-body")
            ], className="card mb-4 shadow")
        ]

    for _id, sample in samples.iterrows():
        try:
            reads_script += "#{}\nln -s {} .\nln -s {} .\n".format(
                sample["name"],
                sample["properties.datafiles.summary.paired_reads"][0],
                sample["properties.datafiles.summary.paired_reads"][1])
        except KeyError as e:
            reads_errors.append("Missing data for sample: {} - {}. In database:\n{}".format(
                sample.get("name", "None"),
                sample["_id"],
                sample.get(
                    "properties.datafiles.summary.paired_reads", "No data")
            ))
    if len(reads_errors):
        reads_html = [
            html.Div([
                html.Div([
                    html.H6("Read files",
                            className="m-0 font-weight-bold text-primary"),
                ], className="card-header py-3"),
                html.Div([
                    "A few errors occurred locating the read paths. If you need more info, " +
                    "please contact an admin.",
                    html.Pre("\n".join(reads_errors), className="error-pre"),
                    html.Pre(reads_script, className="folder-pre")
                ], className="card-body")
            ], className="card shadow")
        ]
    else:
        reads_html = [
            html.Div([
                html.Div([
                    html.H6("Read files",
                            className="m-0 font-weight-bold text-primary"),
                ], className="card-header py-3"),
                html.Div([
                    html.Pre(reads_script, className="folder-pre")
                ], className="card-body")
            ], className="card shadow")
        ]
    return html.Div([html.Div(assemblies_html), html.Div(reads_html)])
