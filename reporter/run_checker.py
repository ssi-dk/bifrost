# -*- coding: utf-8 -*-
import os
import sys
import subprocess

import dash
import dash_auth
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import plotly.figure_factory as ff
from dash.dependencies import Input, Output, State
import keys

import dash_scroll_up

import components.mongo_interface
import components.import_data as import_data


app = dash.Dash()

app.config["suppress_callback_exceptions"] = True
app.title = "bifrost Run Checker"

if hasattr(keys, "pass_protected") and keys.pass_protected:
    dash_auth.BasicAuth(
        app,
        keys.USERNAME_PASSWORD
    )

# Temp css to make it look nice
# Dash CSS
app.css.append_css(
    {"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})
# Lato font
app.css.append_css(
    {"external_url": "https://fonts.googleapis.com/css?family=Lato"})



app.layout = html.Div([
    html.Div(className="container", children=[
        dash_scroll_up.DashScrollUp(
            id="input",
            label="UP",
            className="button button-primary no-print"
        ),
        html.Div(id='none', children=[], style={'display': 'none'}),
        html.H1("SerumQC Run Checker"),
        dcc.Interval(
            id='table-interval',
            interval=30*1000,  # in milliseconds
            n_intervals=0),
        html.Div([
            html.Div(
                [
                    html.Div(
                        dcc.Dropdown(
                            id="run-list",
                            options=[],
                            placeholder="Sequencing run"
                        )
                    )
                ],
                className="nine columns"
            ),
            html.Div(
                [
                    dcc.Link(
                        "Load run",
                        id="run-link",
                        href="/",
                        className="button button-primary")
                ],
                className="three columns"
            )
        ], id="run-selector", className="row"),
        html.P("", id="run-name", style={"display": "none"}),
        html.H2("", id="run-name-title"),
        html.Div(id="report-link"),
        dcc.Location(id="url", refresh=False),
        html.Div(id="run-report", className="run_checker_report"),
        html.Div(id="rerun-form"),
        dcc.ConfirmDialog(message="", displayed=False, id="rerun-output")

    ]),
    html.Footer([
        "Created with 🔬 at SSI. Bacteria icons from ",
        html.A("Flaticon", href="https://www.flaticon.com/"),
        "."], className="footer container")
], className="appcontainer")


# Callbacks

# We could make this one much faster by hiding the unused species with CSS
# by adding a new hidden class.


@app.callback(
    Output("run-name", "children"),
    [Input("url", "pathname")]
)
def update_run_name(pathname):
    if pathname is None:
        pathname = "/"
    path = pathname.split("/")
    if import_data.check_run_name(path[1]):
        name = path[1]
    elif path[1] == "":
        name = ""
    else:
        name = "Not found"
    if len(path) > 2:
        return name + "/" + path[2]
    else:
        return name + "/"


@app.callback(
    Output("run-list", "options"),
    [Input("none", "children")]
)
def update_run_options(none):
    run_list = import_data.get_run_list()
    return [
        {
            "label": "{} ({})".format(run["name"],
                                    len(run["samples"])),
            "value": run["name"]
        } for run in run_list]

@app.callback(
    Output("run-name-title", "children"),
    [Input("run-name", "children")]
)
def update_run_title(run_name):
    return run_name.split("/")[0]


@app.callback(
    Output("report-link", "children"),
    [Input("run-name", "children")]
)
def update_run_name(run_name):
    run_name = run_name.split("/")[0]
    if run_name == "" or run_name == "Not found":
        return None
    else:
        return html.H3(html.A("QC Report", href="{}/{}".format(keys.qc_report_url, run_name)))


@app.callback(
    Output("run-link", "href"),
    [Input("run-list", "value")]
)
def update_run_button(run):
    if run is not None:
        return "/" + run
    else:
        return "/"


@app.callback(
    Output("run-report", "children"),
    [Input("run-name", "children"),
     Input("table-interval", "n_intervals")]
)
def update_run_report(run, n_intervals):
    update_notice = "The table will update every 30s automatically."
    split = run.split("/")
    run = split[0]
    report = "status"
    if len(split) > 1 and split[1] == "resequence":
        report = "resequence"
    data = []
    figure = {}
    if run != "" and report == "status":
        run_data = import_data.get_run(run)
        resequence_link = html.H4(html.A(
            "Resequence Report", href="/{}/resequence".format(run)))
        if run_data is not None:
            components = list(
                map(lambda x: x["name"], run_data["components"]))
        else:
            components = []
        samples = import_data.get_sample_component_status(run)
        header = html.Tr([
            html.Th(html.Div(html.Strong("Sample")),
                    className="rotate rotate-short"),
            html.Th(html.Div(html.Strong("QC status")),
                    className="rotate rotate-short"),
            html.Th()
            ] + list(map(lambda x: html.Th(html.Div(html.Strong(x)), className="rotate rotate-short"), components)))
        rows = [header]
        for name, s_components in samples.items():
            if name == "Undetermined":
                continue
            row = []
            row.append(html.Td(name))
            stamps = import_data.get_sample(
                str(s_components["sample._id"])).get("stamps", {})
            qc_val = stamps.get("ssi_stamper", {}).get("value", "N/A")

            expert_check = False
            if "ssi_expert_check" in stamps and "value" in stamps["ssi_expert_check"]:
                qc_val = stamps["ssi_expert_check"]["value"]
                expert_check = True

            statusname = ""
            if qc_val == "fail:supplying lab":
                qc_val = "SL"
                statusname = "status-1"
            elif qc_val == "N/A":
                statusname = "status--2"
            elif qc_val == "fail:core facility":
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

            for component in components:
                if component in s_components.keys():
                    s_c = s_components[component]
                    row.append(
                        html.Td(s_c[1],
                                className="center status-{}".format(s_c[0])))
                else:
                    row.append(html.Td("None", className="center status-0"))
            rows.append(html.Tr(row))
        table = html.Table(rows, className="unset-width-table")
        update_notice += " Req.: Requirements not met. Init.: initialised."

        return [
            resequence_link,
            html.P(update_notice),
            table]

    if run != "" and report == "resequence":

        run_checker_link = html.H4(html.A(
            "Run Checker Report", href="/{}".format(run)))

        last_runs = import_data.get_last_runs(run, 12)
        last_runs_names = [run["name"] for run in last_runs]
        prev_runs_dict = import_data.get_sample_QC_status(last_runs)
        header = html.Tr([html.Th(html.Div(html.Strong("Sample")), className="rotate")] +
                         list(map(lambda x: html.Th(html.Div(html.Strong(x)), className="rotate"), last_runs_names)))
        rows = [header]
        for name, p_runs in prev_runs_dict.items():
            if name == "Undetermined":
                continue
            row = []
            row.append(html.Td(name))
            
            sample_all_OKs = True

            for index in range(len(last_runs)):
                if last_runs[index]["name"] in p_runs:
                    className = "0"
                    title = "Not Run"
                    status = p_runs[last_runs[index]["name"]]
                    if status.startswith("OK"):
                        className = "2"
                        title = "OK"
                    elif status == "SL":
                        sample_all_OKs = False
                        className = "1"
                        title = "Supplying Lab"
                    elif status.startswith("CF"):
                        # Wont be triggered by supplying because its after.
                        className = "-1"
                        sample_all_OKs = False
                        title = "Core Facility"
                    else:
                        #to account for libray fails
                        sample_all_OKs = False
                    row.append(
                        html.Td(status,
                                className="center status-" + className))
                else:
                    row.append(html.Td("-", className="center status-0"))


            if not sample_all_OKs:
                rows.append(html.Tr(row))
        table = html.Table(rows, className="unset-width-table")
        update_notice += " SL: Supplying Lab, CF: Core Facility, CF(LF): Core Facility (Library Fail). -: No data."
        return [
            run_checker_link,
            html.P(update_notice),
            table
        ]
    return []


@app.callback(
    Output("rerun-form", "children"),
    [Input("run-name", "children")]
)
def update_rerun_form(run_name):
    run_name = run_name.split("/")[0]
    if run_name == "" or not hasattr(keys, "rerun"):
        return None

    run_data = import_data.get_run(run_name)
    if run_data is not None:
        components = list(
            map(lambda x: x["name"], run_data["components"]))
    else:
        components = []

    samples = import_data.get_sample_component_status(run_name)
    samples_options = [{"value": s, "label": s} for s in samples.keys()]
    components_options = [{"value": c, "label": c} for c in components]

    return html.Div([
        html.H6("Rerun components"),
        html.Div([
            html.Div([
                dcc.Dropdown(
                    id="rerun-samples",
                    options=samples_options,
                    placeholder="Samples",
                    multi=True,
                    value=""
                )
            ], className="four columns"),
            html.Div([
                dcc.Dropdown(
                    id="rerun-components",
                    options=components_options,
                    placeholder="Components",
                    multi=True,
                    value=""
                )
            ], className="four columns"),
            html.Div([
                html.Button(
                    "Send",
                    id="rerun-button",
                    n_clicks=0
                )
            ], className="two columns")
        ], className="row"),
    ])

@app.callback(
    Output("rerun-output", "message"),
    [Input("rerun-button", "n_clicks")],
    [State("rerun-samples", "value"),
     State("rerun-components", "value"),
     State("run-name", "children")]
)
def rerun_form_button(button, samples, components, run_name):
    if button == 0 or not hasattr(keys, "rerun") or run_name == "/":
        return ""
    out = []
    run_name = run_name.split("/")[0]
    run_path_bifrost = import_data.get_run(run_name).get("path", "") # ends in /bifrost
    run_path = os.path.dirname(run_path_bifrost) # removes /bifrost
    for sample in samples:
        for component in components:
            command = r'if [ -d \"{}\" ]; then rm -r {}; fi; '.format(
                component, component)
            # unlock first
            command += r"snakemake --shadow-prefix /scratch --restart-times 2 --cores 4 -s {}/src/components/{}/pipeline.smk --config Sample=sample.yaml --unlock; ".format(
                run_path, component)
            command += r"snakemake --shadow-prefix /scratch --restart-times 2 --cores 4 -s {}/src/components/{}/pipeline.smk --config Sample=sample.yaml".format(
                run_path, component)
            if keys.rerun["grid"] == "slurm":
                process = subprocess.Popen('sbatch --mem={memory}G -p {priority} -c {threads} -t {walltime} -J "bifrost_{sample_name}" --wrap "{command}"'.format(
                    **keys.rerun, sample_name=sample, command=command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=os.environ, cwd=run_path + "/" + sample)
                process_out, process_err = process.communicate()
                out.append((sample, component, process_out, process_err))
            elif keys.rerun["grid"] == "torque":
                
                if "advres" in keys.rerun:
                    advres = ",advres={}".format(
                        keys.rerun["advres"])
                else:
                    advres = ''
                torque_node = ",nodes=1:ppn={}".format(keys.rerun["threads"])
                script_path = os.path.join(run_path, sample, "manual_rerun.sh")
                with open(script_path, "w") as script:
                    command += "#PBS -V -d . -w . -l mem={memory}gb,nodes=1:ppn={threads},walltime={walltime}{advres} -N 'bifrost_{sample_name}' -W group_list={group} -A {group} \n".format(
                        **keys.rerun, sample_name=sample
                    )
                    script.write(command)
                process = subprocess.Popen('qsub {}'.format(script_path), stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT, shell=True, env=os.environ, cwd=run_path + "/" + sample)
                process_out, process_err = process.communicate()
                out.append((sample, component, process_out, process_err))

    message = "Jobs sent to the server:\n" + "\n".join(["{}, {}: out: {} | err: {}".format(*el) for el in out])
    message += "\nClick OK or Cancel to close this notice."
    return message


@app.callback(
    Output("rerun-output", "displayed"),
    [Input("rerun-output", "message")],
)
def update_run_report(message):
    if message != "":
        return True
    return False

server = app.server  # Required for gunicorn

#app.run_server(debug=True, host="0.0.0.0")
if __name__ == '__main__':
    # 0.0.0.0 exposes the app to the network.
    app.run_server(debug=True, host="0.0.0.0",
                   port=8051, dev_tools_hot_reload=True)
