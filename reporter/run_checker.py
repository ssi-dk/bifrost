# -*- coding: utf-8 -*-
import os
import sys
import subprocess

import dash
import dash_auth
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objs as go
import plotly.figure_factory as ff
from dash.dependencies import Input, Output, State
import keys

import dash_scroll_up

import components.mongo_interface
import components.import_data as import_data

import json
from bson import json_util



def pipeline_report(sample_data):
    update_notice = "The table will update every 30s automatically."


    table = dash_table.DataTable(
        id="pipeline-table", selected_rows=[],
        style_table={
            'overflowX': 'scroll',
        },
        page_action='none',
        style_as_list_view=True,
        style_cell={
            'textAlign': 'center',
            "fontFamily": "Arial",
            "padding": "0px 10px",
            "fontSize": "0.7rem",
            "height": "25px"
            })

    update_notice += (" Req.: requirements not met. Init.: initialised. "
                      "*: user submitted")
    rerun_columns = [
        {"id": "sample", "name": "sample"},
        {"id": "component", "name": "component"},
    ]

    return [
        # resequence_link,
        dbc.Row([
            dbc.Col([
                html.H2("Pipeline Status", className="mt-3"),
                html.P(update_notice),
                table,
                dcc.Interval(
                    id='table-interval',
                    interval=30*1000,  # in milliseconds
                    n_intervals=0
                )
            ], width=9),
            dbc.Col([
                html.H3("Rerun components", className="mt-3"),
                dbc.Alert(id="rerun-output",
                          color="secondary",
                          dismissable=True,
                          is_open=False),
                html.Label(html.Strong("Add all components for sample")),
                dbc.InputGroup([
                    dcc.Dropdown(id="rerun-components",
                                className="dropdown-group"),
                    dbc.InputGroupAddon(
                        dbc.Button("Add", id="rerun-add-components"),
                        addon_type="append",
                    ),
                ]),
                html.Label(html.Strong("Add component for all samples"),
                           className="mt-3"),
                dbc.InputGroup([
                    dcc.Dropdown(id="rerun-samples",
                                 className="dropdown-group"),
                    dbc.InputGroupAddon(
                        dbc.Button("Add", id="rerun-add-samples"),
                        addon_type="append",
                    ),
                ]),
                html.Label(html.Strong("Add all failed components"),
                           className="mt-3"),
                dbc.Button("Add", id="rerun-add-failed", block=True),
                html.Div([
                    html.H4("Selected sample components"),
                    dash_table.DataTable(
                        id="pipeline-rerun",
                        columns=rerun_columns,
                        row_deletable=True),
                    dbc.Button("Rerun selected sample components",
                               id="rerun-button", className="mt-3", block=True,
                               color="primary"),
                ], className="mt-3")
            ],
                width=3,
                style={"backgroundColor": "rgba(0, 0, 0, .05)"}
            )
        ])
    ]


def pipeline_report_data(sample_data):

    s_c_status = import_data.get_sample_component_status(
        sample_data)

    components_order = [
        "whats_my_species", "analyzer", "assemblatron", "ssi_stamper",
        "ariba_resfinder", "ariba_mlst", "ariba_plasmidfinder",
        "ariba_virulencefinder", "sp_cdiff_fbi", "sp_ecoli_fbi",
        "sp_salm_fbi", "min_read_check", "qcquickie"]


    s_c_components = []
    for sample_id, s_c in s_c_status.items():
        for comp_name in s_c.keys():
            if comp_name not in s_c_components and comp_name != "sample":
                s_c_components.append(comp_name)
    components = [comp for comp in components_order if comp in s_c_components]

    rows = []

    columns = [
        {"name": "Priority", "id": "priority"},
        {"name": "Sample", "id": "sample"},
        {"name": "QC status", "id": "qc_val"}
    ]

    rerun_form_components = []

    for comp in components:
        columns.append({"name": comp, "id": comp})
        rerun_form_components.append({"label": comp, "value": comp})

    # Conditional data colors
    style_data_conditional = [
        {
            "if": {
                "column_id": "qc_val",
                "filter": 'qc_val eq "CF"'
            },
            "backgroundColor": "#ea6153"
        },
        {
            "if": {
                "column_id": "qc_val",
                "filter": 'qc_val eq "CF(LF)"'
            },
            "backgroundColor": "#ea6153"
        },
        {
            "if": {
                "column_id": "qc_val",
                "filter": 'qc_val eq "OK"'
            },
            "backgroundColor": "#27ae60"
        },
        {
            "if": {
                "column_id": "qc_val",
                "filter": 'qc_val eq "SL"'
            },
            "backgroundColor": "#f1c40f"
        }
    ]
    for col in s_c_components:
        style_data_conditional.append({
            "if": {
                "column_id": col,
                "filter": '{} eq "Fail"'.format(col)
            },
            "backgroundColor": "#ea6153"
        })
        style_data_conditional.append({
            "if": {
                "column_id": col,
                "filter": '{} eq "OK"'.format(col)
            },
            "backgroundColor": "#3498db"
        })
        style_data_conditional.append({
            "if": {
                "column_id": col,
                "filter": '{} eq "Running"'.format(col)
            },
            "backgroundColor": "#f1c40f"
        })
        style_data_conditional.append({
            "if": {
                "column_id": col,
                "filter": '{} eq "Req."'.format(col)
            },
            "backgroundColor": "#d3d3d3",
            "color": "#525252"
        })

    rerun_form_samples = []

    for sample_id, s_components in s_c_status.items():
        row = {}
        sample = s_components["sample"]

        name = sample["name"]
        row["sample"] = name
        row["_id"] = str(sample["_id"])
        if name == "Undetermined":
            continue  # ignore this row

        rerun_form_samples.append({"label": name, "value": "{}:{}".format(
            sample_id, name)})

        stamps = sample.get("stamps", {})
        priority = sample.get("sample_sheet",
                              {}).get("priority", "").lower()
        prio_display = " "
        prio_title = priority
        if priority == "high":
            prio_display = "🚨"
        else:
            prio_display = ""
        row["priority"] = prio_display
        qc_val = stamps.get("ssi_stamper", {}).get("value", "N/A")

        expert_check = False
        if ("supplying_lab_check" in stamps and
                "value" in stamps["supplying_lab_check"]):
            qc_val = stamps["supplying_lab_check"]["value"]
            expert_check = True

        statusname = ""
        if qc_val == "fail:supplying lab":
            qc_val = "SL"
            statusname = "status-1"
        elif qc_val == "N/A":
            statusname = "status--2"
        elif (qc_val == "fail:core facility" or
                qc_val == "fail:resequence"):
            statusname = "status--1"
            qc_val = "CF"
        elif qc_val == "pass:OK":
            statusname = "status-2"
            qc_val = "OK"

        if expert_check:
            qc_val += "*"

        row["qc_val"] = qc_val

        for component in components:
            if component in s_components.keys():
                s_c = s_components[component]
                row[component] = s_c[1]
            else:
                row[component] = "None"
        rows.append(row)
    return rows, columns, style_data_conditional, rerun_form_components, rerun_form_samples

# # Callbacks

# @app.callback(
#     Output("run-name", "children"),
#     [Input("url", "pathname")]
# )
# def update_run_name(pathname):
#     if pathname is None:
#         pathname = "/"
#     path = pathname.split("/")
#     if import_data.check_run_name(path[1]):
#         name = path[1]
#     elif path[1] == "":
#         name = ""
#     else:
#         name = "Not found"
#     if len(path) > 2:
#         return name + "/" + path[2]
#     else:
#         return name + "/"


# @app.callback(
#     Output("run-list", "options"),
#     [Input("none", "children")]
# )
# def update_run_options(none):
#     run_list = import_data.get_run_list()
#     return [
#         {
#             "label": "{} ({})".format(run["name"],
#                                       len(run["samples"])),
#             "value": run["name"]
#         } for run in run_list]


# @app.callback(
#     Output("run-name-title", "children"),
#     [Input("run-name", "children")]
# )
# def update_run_title(run_name):
#     return run_name.split("/")[0]


# @app.callback(
#     Output("report-link", "children"),
#     [Input("run-name", "children")]
# )
# def update_run_name(run_name):
#     run_name = run_name.split("/")[0]
#     if run_name == "" or run_name == "Not found":
#         return None
#     else:
#         return html.H3(html.A("QC Report",
#                               href="{}/{}".format(keys.qc_report_url,
#                                                   run_name)))


# @app.callback(
#     Output("run-link", "href"),
#     [Input("run-list", "value")]
# )
# def update_run_button(run):
#     if run is not None:
#         return "/" + run
#     else:
#         return "/"


# @app.callback(
#     Output("sample-store", "data"),
#     [Input("run-name", "children")]
# )
# def update_sample_store(run_name):
#     split = run_name.split("/")
#     run_name = split[0]
#     run = import_data.get_run(run_name)
#     if run is None:
#         return json_util.dumps({"run": None})
#     samples = import_data.get_samples(list(map(lambda x: str(x["_id"]), run["samples"])))

#     store = {
#         "run": run,
#         "samples_by_id": {str(s["_id"]): s for s in samples},
#         "report": "resequence" if len(split) > 1 and split[1] == "resequence" else "run_checker"
#     }

#     return json_util.dumps(store)


# @app.callback(
#     Output("run-report", "children"),
#     [Input("sample-store", "data"),
#      Input("table-interval", "n_intervals")]
# )
# def update_run_report(store, n_intervals):
#     update_notice = "The table will update every 30s automatically."
#     store = json_util.loads(store)
#     run = store["run"]
#     if run is None:
#         return None
#     samples_by_id = store["samples_by_id"]

#     if store["report"] == "run_checker":
#         run_data = import_data.get_run(run)  #NOTE: de;ete tjos
        
#         resequence_link = html.H4(dcc.Link(
#             "Resequence Report", href="/{}/resequence".format(run["name"]))) #move to multiple outputs

#         s_c_status = import_data.get_sample_component_status(samples_by_id.keys())

#         if "components" in run:
#             components = list(
#                 map(lambda x: x["name"], run["components"]))
#         else:
#             components = []
#             for sample_id, s_c in s_c_status.items():
#                 if s_c["component"]["name"] not in components:
#                     components.append(s_c["component"]["name"])
#         header = html.Tr([
#             html.Th(html.Div(html.Strong("Priority")),
#                     className="rotate rotate-short"),
#             html.Th(html.Div(html.Strong("Sample")),
#                     className="rotate rotate-short"),
#             html.Th(html.Div(html.Strong("QC status")),
#                     className="rotate rotate-short"),
#             html.Th()
#             ] + list(map(lambda x: html.Th(html.Div(html.Strong(x)),
#                                            className="rotate rotate-short"),
#                          components)))
#         rows = [header]
#         for sample_id, s_components in s_c_status.items():
#             sample = samples_by_id[sample_id]
#             name = sample["name"]
#             if name == "Undetermined":
#                 continue #ignore this row
#             row = []
#             sample = import_data.get_sample(
#                 str(s_components["sample._id"]))
#             stamps = sample.get("stamps", {})
#             priority = sample.get("sample_sheet",
#                                   {}).get("priority", "").lower()
#             prio_display = " "
#             prio_title = priority
#             if priority == "high":
#                 prio_display = "🚨"
#             else:
#                 prio_display = ""
#             row.append(html.Td(prio_display,
#                                title=prio_title,
#                                className="center"))

#             row.append(html.Td(name))
#             qc_val = stamps.get("ssi_stamper", {}).get("value", "N/A")

#             expert_check = False
#             if ("supplying_lab_check" in stamps and
#                     "value" in stamps["supplying_lab_check"]):
#                 qc_val = stamps["supplying_lab_check"]["value"]
#                 expert_check = True

#             statusname = ""
#             if qc_val == "fail:supplying lab":
#                 qc_val = "SL"
#                 statusname = "status-1"
#             elif qc_val == "N/A":
#                 statusname = "status--2"
#             elif (qc_val == "fail:core facility" or
#                   qc_val == "fail:resequence"):
#                 statusname = "status--1"
#                 qc_val = "CF"
#             elif qc_val == "pass:OK":
#                 statusname = "status-2"
#                 qc_val = "OK"

#             if expert_check:
#                 qc_val += "*"

#             row.append(
#                 html.Td(qc_val, className="center {}".format(statusname)))
#             row.append(html.Td())

#             for component in components:
#                 if component in s_components.keys():
#                     s_c = s_components[component]
#                     row.append(
#                         html.Td(s_c[1],
#                                 className="center status-{}".format(s_c[0])))
#                 else:
#                     row.append(html.Td("None", className="center status-0"))
#             rows.append(html.Tr(row))
#         table = html.Table(rows, className="unset-width-table")
#         update_notice += (" Req.: requirements not met. Init.: initialised. "
#                           "*: user submitted")

#         return [
#             resequence_link,
#             html.P(update_notice),
#             table]

#     if store["report"] == "resequence":

#         run_checker_link = html.H4(html.A(
#             "Run Checker Report", href="/{}".format(run)))

#         last_runs = import_data.get_last_runs(run["name"], 12) #Get last 12 runs
#         last_runs_names = [run["name"] for run in last_runs]
#         prev_runs_dict = import_data.get_sample_QC_status(last_runs)
#         header = html.Tr([html.Th(html.Div(html.Strong("Sample")),
#                                   className="rotate")] +
#                          list(map(lambda x: html.Th(html.Div(html.Strong(x)),
#                                                     className="rotate"),
#                                   last_runs_names)))
#         rows = [header]
#         for name, p_runs in prev_runs_dict.items():
#             if name == "Undetermined":
#                 continue
#             row = []
#             row.append(html.Td(name))

#             sample_all_OKs = True

#             for index in range(len(last_runs)):
#                 if last_runs[index]["name"] in p_runs:
#                     className = "0"
#                     title = "Not Run"
#                     status = p_runs[last_runs[index]["name"]]
#                     if status.startswith("OK"):
#                         className = "2"
#                         title = "OK"
#                     elif status == "SL":
#                         sample_all_OKs = False
#                         className = "1"
#                         title = "Supplying Lab"
#                     elif status.startswith("CF"):
#                         # Wont be triggered by supplying because its after.
#                         className = "-1"
#                         sample_all_OKs = False
#                         title = "Core Facility"
#                     else:
#                         # to account for libray fails
#                         sample_all_OKs = False
#                     row.append(
#                         html.Td(status,
#                                 className="center status-" + className))
#                 else:
#                     row.append(html.Td("-", className="center status-0"))

#             if not sample_all_OKs:
#                 rows.append(html.Tr(row))
#         table = html.Table(rows, className="unset-width-table")
#         update_notice += (" SL: Supplying Lab, CF: Core Facility, CF(LF): "
#                           "Core Facility (Library Fail). -: No data. "
#                           "*: user submitted")
#         return [
#             run_checker_link,
#             html.P(update_notice),
#             table
#         ]
#     return []


# @app.callback(
#     Output("rerun-form", "children"),
#     [Input("run-name", "children")]
# )
# def update_rerun_form(run_name):
#     run_name = run_name.split("/")[0]
#     if run_name == "" or not hasattr(keys, "rerun"):
#         return None

#     run_data = import_data.get_run(run_name)
#     if run_data is not None:
#         components = list(
#             map(lambda x: x["name"], run_data["components"]))
#     else:
#         components = []

#     components_options = [{"value": c, "label": c} for c in components]

#     return html.Div([
#         html.H6("Rerun components"),
#         html.Div([
#             html.Div([
#                 dcc.Textarea(
#                     id="rerun-samples",
#                     placeholder="one sample per line",
#                     value=""
#                 )
#             ], className="four columns"),
#             html.Div([
#                 dcc.Dropdown(
#                     id="rerun-components",
#                     options=components_options,
#                     placeholder="Components",
#                     multi=True,
#                     value=""
#                 )
#             ], className="four columns"),
#             html.Div([
#                 html.Button(
#                     "Send",
#                     id="rerun-button",
#                     n_clicks=0
#                 )
#             ], className="two columns")
#         ], className="row"),
#     ])



def rerun_components_button(button, table_data):
    if button == 0:
        return "", False
    out = []
    to_rerun = {}
    for row in table_data:
        sample_rerun = to_rerun.get(row["sample_id"], [])
        sample_rerun.append(row["component"])
        to_rerun[row["sample_id"]] = sample_rerun
    
    sample_dbs = import_data.get_samples(sample_ids=to_rerun.keys())
    samples_by_id = {str(s["_id"]) : s for s in sample_dbs}

    bifrost_components_dir = os.path.join(keys.rerun["bifrost_dir"], "components/")

    for sample, components in to_rerun.items():
        sample_db = samples_by_id[sample]
        sample_name = sample_db["name"]
        run_path = sample_db["path"]
        sample_command = ""
        for component in components:
            component_path = os.path.join(bifrost_components_dir,
                                          component, "pipeline.smk")
            command = r'if [ -d \"{}\" ]; then rm -r {}; fi; '.format(
                component, component)
            # unlock first
            command += (r"snakemake --shadow-prefix /scratch --restart-times 2"
                        r" --cores 4 -s {} "
                        r"--config Sample=sample.yaml --unlock; ").format(
                            component_path)
            command += (r"snakemake --shadow-prefix /scratch --restart-times 2"
                        r" --cores 4 -s {} "
                        r"--config Sample=sample.yaml; ").format(
                            component_path)
            sample_command += command
        
        if keys.rerun["grid"] == "slurm":
            process = subprocess.Popen(
                ('sbatch --mem={memory}G -p {priority} -c {threads} '
                    '-t {walltime} -J "bifrost_{sample_name}" --wrap'
                    ' "{command}"').format(
                        **keys.rerun,
                        sample_name=sample_name,
                        command=sample_command),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                env=os.environ,
                cwd=run_path)
            process_out, process_err = process.communicate()
            out.append((sample_name, process_out, process_err))
        elif keys.rerun["grid"] == "torque":

            if "advres" in keys.rerun:
                advres = ",advres={}".format(
                    keys.rerun["advres"])
            else:
                advres = ''
            torque_node = ",nodes=1:ppn={}".format(keys.rerun["threads"])
            script_path = os.path.join(run_path, "manual_rerun.sh")
            with open(script_path, "w") as script:
                command += ("#PBS -V -d . -w . -l mem={memory}gb,nodes=1:"
                            "ppn={threads},walltime={walltime}{advres} -N "
                            "'bifrost_{sample_name}' -W group_list={group}"
                            " -A {group} \n").format(**keys.rerun,
                                                     sample_name=sample_name)
                script.write(command)
            process = subprocess.Popen('qsub {}'.format(script_path),
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT,
                                        shell=True,
                                        env=os.environ,
                                        cwd=run_path)
            process_out, process_err = process.communicate()
            out.append((sample_name, process_out, process_err))
        elif keys.rerun["grid"] == "slurm.mock":
            print(('sbatch --mem={memory}G -p {priority} -c {threads} '
                  '-t {walltime} -J "bifrost_{sample_name}" --wrap'
                  ' "{command}"').format(
                **keys.rerun,
                sample_name=sample_name,
                command=sample_command))

    message = "Jobs sent to the server:\n"
    message += "\n".join(["{}: out: {} | err: {}".format(*el)
                         for el in out])
    return message, True


def update_rerun_table(active, table_data, n_click_comp, n_click_samp,
                       n_click_fail, columns, prev_data, rerun_comp,
                       rerun_samp):
    # default values
    if prev_data is None:
        prev_data = []

    if columns:
        columns = columns[3:]

    #Get context to know which button was triggered.
    ctx = dash.callback_context

    if not ctx.triggered:
        triggered_id = None
    else:
        triggered_id = ctx.triggered[0]['prop_id']

    #Nothing triggered it, return empty table if init call or prev data.
    if active is None and triggered_id is None:
        return prev_data

    if triggered_id == "pipeline-table.active_cell":
        col_index = active[1] - 3
        col = columns[col_index]["id"]
        sample = table_data[active[0]]["sample"]
        sample_id = table_data[active[0]]["_id"]

        new_rows = [{"sample": sample, "component": col,
                     "sample_id": sample_id}]
    elif triggered_id == "rerun-add-components.n_clicks":
        sample_id, sample = rerun_comp.split(":")
        new_rows = [{"sample": sample,
                     "component": comp["id"],
                     "sample_id": sample_id} for comp in columns]
    elif triggered_id == "rerun-add-samples.n_clicks":
        new_rows = []
        for row in table_data:
            new_rows.append({"sample": row["sample"],
                             "component": rerun_samp,
                             "sample_id": row["_id"]})
    elif triggered_id == "rerun-add-failed.n_clicks":
        new_rows = []
        for row in table_data:
            for col in columns:
                col = col["id"]
                if row[col] == "Fail":
                    new_rows.append({"sample": row["sample"],
                                     "component": col,
                                     "sample_id": row["_id"]})
    else:
        new_rows = []

    for new_row in new_rows:
        if new_row not in prev_data:
            prev_data = prev_data + [new_row]

    return prev_data
