# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import urllib.parse as urlparse

import dash
import dash_auth
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import dash_bootstrap_components as dbc
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import plotly.figure_factory as ff
from dash.dependencies import Input, Output, State
import keys

import dash_scroll_up

import components.import_data as import_data

import json
from bson import json_util



def pipeline_report(sample_data):
    update_notice = "The table will update every 30s automatically."

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
                html.Div([
                    html.Div([
                        html.H6("Pipeline status", className="m-0 font-weight-bold text-primary")
                    ], className="card-header py-3"),
                    html.Div([
                        html.P(update_notice),
                        html.Div(id="pipeline-table"),
                        dcc.Interval(
                            id='table-interval',
                            interval=30*1000,  # in milliseconds
                            n_intervals=0
                        )
                    ], className="card-body")
                ], className="card shadow mb-4")
            ], width=9),
            dbc.Col([
                html.Div([
                    html.Div([
                        html.H6("Rerun components",
                                className="m-0 font-weight-bold text-primary")
                    ], className="card-header py-3"),
                    html.Div([
                        dbc.Alert(id="rerun-output",
                                  color="secondary",
                                  dismissable=True,
                                  is_open=False),
                        html.Label(html.Strong(
                            "Add all components for sample")),
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
                    ], className="card-body")
                ], className="card shadow mb-4")
            ], width=3
            )
        ])
    ]


def pipeline_report_data(sample_data):

    status_dict = {
        "Success": "OK",
        "Running": "Running",
        "initialized": "init.",
        "Failure": "Fail",
        "Requirements not met": "Req.",
        "queued to run": "queue",
    }
    status_code_dict = {
        "Success": 2,
        "Running": 1,
        "initialized": 0,
        "Failure": -1,
        "Requirements not met": -2,
        "queued to run": 0
    }
    samples = import_data.filter_all(
        sample_ids = [s["_id"] for s in sample_data],
        include_s_c=False,
        projection={"stamps": 1,
         'reads.R1': 1,
         'name': 1,
         'sample_sheet.priority': 1})
    samples["_id"] = samples["_id"].astype(str)
    samples = samples.set_index("_id")

    s_c_status = import_data.get_sample_component_status(
        sample_data)

    components_order = [
        "whats_my_species", "analyzer", "assemblatron", "ssi_stamper",
        "ariba_resfinder", "ariba_mlst", "ariba_plasmidfinder",
        "ariba_virulencefinder", "sp_cdiff_fbi", "sp_ecoli_fbi",
        "sp_salm_fbi", "min_read_check", "qcquickie"]


    s_c_components = []
    for sample in s_c_status:
        for comp_name in sample['s_cs'].keys():
            if comp_name not in s_c_components:
                s_c_components.append(comp_name)
    components_list = [comp for comp in components_order if comp in s_c_components]

    header = html.Tr([
        html.Th(html.Div(html.Strong("Priority")),
                className="rotate rotate-short"),
        html.Th(html.Div(html.Strong("Sample")),
                className="rotate rotate-short"),
        html.Th(html.Div(html.Strong("QC status")),
                className="rotate rotate-short"),
        html.Th()
    ] + list(map(lambda x: html.Th(html.Div(html.Strong(x)),
                                   className="rotate rotate-short"),
                 components_list)))
    rows = [header]

    rerun_form_components = []

    for comp in components_list:
        rerun_form_components.append({"label": comp, "value": comp})

    rerun_form_samples = []

    for s_components in s_c_status:
        sample_id = str(s_components["_id"])
        sample = samples.loc[sample_id]
        name = sample["name"]

        row = []
        if name == "Undetermined":
            continue  # ignore this row


        rerun_form_samples.append({"label": name, "value": "{}:{}".format(
            sample_id, name)})

        priority = str(sample.get("sample_sheet.priority", "")).lower()
        if priority == "high":
            prio_display = "🚨"
        else:
            prio_display = ""
        prio_title = priority
        row.append(html.Td(prio_display,
                            title=prio_title,
                            className="center"))
        row.append(html.Td(name))
        qc_val = sample.get("stamps.ssi_stamper.value", "N/A")
        if pd.isnull(qc_val):
            qc_val = "N/A"

        if name == '131964':
            print(qc_val)
            print(sample)

        if qc_val == "N/A" and (not pd.isnull(sample.get("reads.R1", np.NaN))):
            qc_val = "CF(LF)"

        expert_check = False
        if sample.get('stamps.supplying_lab_check.value') is not None:
            qc_val = sample.get('stamps.supplying_lab_check.value')
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

        row.append(
            html.Td(qc_val, className="center {}".format(statusname)))
        row.append(html.Td())

        for component in components_list:
            if component in s_components["s_cs"]:
                s_c = s_components["s_cs"][component]
                row.append(
                    html.Td(status_dict[s_c],
                            className="center status-{}".format(status_code_dict[s_c])))
            else:
                row.append(html.Td("None", className="center status-0"))
        rows.append(html.Tr(row))
    table = html.Table(rows, className="unset-width-table")
    return table, rerun_form_components, rerun_form_samples





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
