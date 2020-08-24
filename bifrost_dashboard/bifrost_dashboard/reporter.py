# -*- coding: utf-8 -*-
import os
import re
import urllib.parse as urlparse
import datetime
from flask_caching import Cache

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_auth
import numpy as np
import dash_bootstrap_components as dbc
from dash.dependencies import ClientsideFunction, Input, Output, State
import dash_table


import bifrostapi

import bifrost_dashboard.components.import_data as import_data
from bifrost_dashboard.components.filter import html_div_filter, generate_table, filter_update_run_options, filter_update_filter_values, html_filter_drawer, html_collection_selector, update_collection_button
from bifrost_dashboard.components.sample_report import SAMPLE_PAGESIZE, sample_report, children_sample_list_report, samples_next_page
import bifrost_dashboard.components.global_vars as global_vars
import bifrost_dashboard.components.admin as admin
from bifrost_dashboard.run_checker import pipeline_report, rerun_components_button, update_rerun_table, pipeline_report_data
from bifrost_dashboard.components.aggregate_report import aggregate_report, update_aggregate_fig, aggregate_species_dropdown
from bifrost_dashboard.components.resequence_report import resequence_report
from bifrost_dashboard.components.link_to_files import link_to_files, link_to_files_div
import bifrost_dashboard.components.virtual_runs as virtual_runs

import yaml
config = yaml.safe_load(open(os.environ["BIFROST_DASH_CONFIG"]))


bifrostapi.add_URI(config["mongodb_key"])

external_scripts = [
    'https://kit.fontawesome.com/24170a81ff.js',
]

external_stylesheets = [
    "https://fonts.googleapis.com/css?family=Lato",
    dbc.themes.BOOTSTRAP
]
assets = os.path.dirname(os.path.abspath(__file__)) + "/data/assets"
app = dash.Dash("bifrost_dashboard",
    assets_folder=assets,
    external_stylesheets=external_stylesheets,
    external_scripts=external_scripts
)
app.title = "bifrost"
app.config["suppress_callback_exceptions"] = True
cache = Cache(app.server, config={
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': config["cache_location"]
})
cache_timeout = 60

if config.get("pass_protected"):
    dash_auth.BasicAuth(
        app,
        config.USERNAME_PASSWORD
    )

# Temp css to make it look nice
# Lato font


def samples_list(active, collection_name=None):
    links = [
        {
            "icon": "fa-list",
            "href": ""
        },
        {
            "icon": "fa-money-check",
            "href": "sample-report"
        },
        {
            "icon": "fa-chart-pie",
            "href": "aggregate"
        },
        {
            "icon": "fa-traffic-light",
            "href": "pipeline-report"
        },
        {
            "icon": "fa-link",
            "href": "link-to-files"
        }
    ]
    link_list = []
    for item in links:
        href = "/" + item["href"]
        if collection_name is not None:
            href = "/collection/{}/{}".format(collection_name, item["href"])
        if active == item['href']:
            link_list.append(dcc.Link(
                html.I(className="fas {} fa-fw".format(item['icon'])),
                className="btn btn-outline-secondary active",
                href=href
            ))
        else:
            link_list.append(dcc.Link(
                html.I(className="fas {} fa-fw".format(item['icon'])),
                className="btn btn-outline-secondary",
                href=href
            ))
    return link_list

app.layout = html.Div([
    dcc.Location(id="url", refresh=False),
    # To store url param values
    dcc.Store(id="sample-store", data=[], storage_type='session'),
    dcc.Store(id="param-store", data={}),
    dcc.Store(id="selected-collection", data=None),
    html.Ul(
        [
            html.A(
                [
                    html.Img(src="/assets/img/bifrost-logo-white@2x.png",
                             className="navbar-logo ")
                    # html.Div("bifrost", className="sidebar-brand-text mx-3")
                ],
                className="sidebar-brand d-flex align-items-center justify-content-center",
                href="/"
            ),
            
            html.Hr(className="sidebar-divider"),
            html.Div("Browse data", className="sidebar-heading"),

            html.Li(dcc.Link(
                [
                    html.I(className="fas fa-vial fa-fw"),
                    html.Span("Samples")
                ], className="nav-link", href="/"),
                className="nav-item",
                id="samples-nav"),
            html.Li(dcc.Link([
                    html.I(className="fas fa-vials fa-fw"),
                    html.Span("Collections")
                    ], className="nav-link", href="/collection"),
                    className="nav-item", id="collections-nav"),
            html.Hr(className="sidebar-divider"),
            html.Div("Reports", className="sidebar-heading"),

            html.Li(dcc.Link(
                [
                    html.I(className="fas fa-chart-line fa-fw"),
                    html.Span("Resequence report")
                ], className="nav-link", href="/resequence-report"),
                className="nav-item", id="resequence-nav"),
            html.Hr(className="sidebar-divider"),
            html.Div(
                html.Button(className="rounded-circle border-0",
                            id="sidebarToggle"),
                className="text-center d-none d-md-inline"
            ),
        ],
        className="navbar-nav bg-gradient-primary sidebar sidebar-dark accordion",
        id="sidebar"
    ),
    html.Div([
        html.Div(id="content", children=[
            html.Nav(
                [
                    html.Ul([
                        html.Li(
                            html.Span("This view is in Beta. Please report any feedback/bugs to mbas@ssi.dk :)"), className="nav-item mx-1"
                        ),
                    ], className="navbar-nav"),
                    html.Ul([
                        html.Li(
                            html.A("Documentation", href="https://ssi-dk.github.io/bifrost/"), className="nav-item dropdown no-arrow mx-1"
                        ),
                        html.Li(
                            html.A("Github", href="https://github.com/ssi-dk/bifrost"), className="nav-item dropdown no-arrow mx-1"
                        )
                    ],className="navbar-nav ml-auto")
                ],
                className="navbar navbar-expand navbar-light bg-white topbar mb-4 static-top shadow"),
            html.Main([
                html_collection_selector(),
                html.Div([
                    dbc.Collapse(
                        [
                            html_filter_drawer()
                        ], id="filter_panel"
                    ),
                    html.Div([
                        html.Div([
                            html.Div(
                                samples_list('/'),
                                className="btn-group shadow-sm",
                                id="selected-view-buttons"
                            ),
                        ], className="col-4"),
                        html.Div([
                            html.Button(
                                html.I(className="fas fa-filter fa-sm"),
                                className="btn btn-outline-secondary shadow-sm mx-auto d-block",
                                id="filter_toggle"
                            ),
                        ], className="col-4"),
                    ], className="row mb-4"),
                ], id="samples-panel", className="d-none"),
                html.Div(id="selected-view"),
                virtual_runs.modal_add(),
                virtual_runs.modal_remove(),
                html.Footer([
                    "Created with 🔬 at SSI. Bacteria icons from ",
                    html.A("Flaticon", href="https://www.flaticon.com/"),
                    "."], className="footer container")
            ], className="container-fluid",
            role="main"),
        ]),
    ], id="content-wrapper", className="d-flex flex-column")
], id="wrapper")


# Callbacks

# We could make this one much faster by hiding the unused species with CSS
# by adding a new hidden class.



@app.callback(
    [Output("filter_panel", "is_open"),
     Output("filter_toggle", "className")],
    [Input("filter_toggle", "n_clicks")],
    [State("filter_panel", "is_open")]
)
def sidebar_toggle(n_clicks, is_open):
    if n_clicks:
        if is_open:
            return [False, "btn btn-outline-secondary shadow-sm mx-auto d-block"]
        else:
            return [True, "btn btn-outline-secondary shadow-sm mx-auto d-block active"]
    return [is_open, "btn btn-outline-secondary shadow-sm mx-auto d-block"]


@app.callback(
    Output("param-store", "data"),
    [Input("url", "search")],
    [State("param-store", "data")]
)
def update_run_name(params, prev_params):
    if params is None or params == "":
        raise dash.exceptions.PreventUpdate("Initial repeated params call")
    pparse = urlparse.urlparse(params)
    params = urlparse.parse_qs(pparse.query)
    if params == prev_params:
        raise dash.exceptions.PreventUpdate("No param change")
    return params

@app.callback(
    [Output("selected-view", "children"),
     Output("selected-view-buttons", "children"),
     Output("samples-panel", "className"),
     Output("samples-nav", "className"),
     Output("resequence-nav", "className"),
     Output("collections-nav", "className"),
     Output("selected-collection", "data"),
     Output("collection-selector-div", "className"),
     Output("run-list-div", "className")],
    [Input("url", "pathname")],
    [State("sample-store", "data")]
)
def update_view(pathname, sample_store):

    if pathname is None or pathname == "/":
        pathname = "/"
    path = pathname.split("/")
    view = None
    samples_panel = ""
    samples_nav = "nav-item"
    resequence_nav = "nav-item"
    collections_nav = "nav-item"
    collection_view = False
    collection_name = None
    if path[1] == "collection":
        collection_view = True
        if len(path) > 2: #/collection/collectionname
            collection_name = path[2]
            if len(path) > 3:  # /collection/collectionname/section
                section = path[3]
            else:  # /collection/collectionname
                section = ""
        else:  # /collection
            section = ""
    else:  # /section
        section = path[1]
        if section == "resequence-report":
            if len(path) == 3:  #/resequence-report/collectionname
                collection_name = path[2]
            resequence_nav += " active" 
        else:
            samples_nav += " active"


    if section == "":
        run = bifrostapi.get_run(collection_name)
        if collection_view and run is not None and run["type"] == "virtual":
            view = html_div_filter(True)
        else:
            view = html_div_filter(False)
    elif section == "sample-report":
        view = sample_report(sample_store)
    elif section == "aggregate":
        view = aggregate_report()  # Doesn't need data
    elif section == "pipeline-report":
        view = pipeline_report(sample_store)
    elif section == "resequence-report":
        samples_panel = "d-none"
        view = resequence_report(collection_name)
    elif section == "link-to-files":
        view = link_to_files_div()
    else:
        samples_panel = "d-none"
        view = "Not found"
    
    if collection_view:
        collection_selector_list = "row"
        run_list = "d-none"
        collections_nav += " active"
    elif section == "resequence-report":
        collection_selector_list = "row"
        run_list = "d-none"
    else:
        collection_selector_list = "row d-none"
        run_list = ""

    return [view, samples_list(section, collection_name), samples_panel,
            samples_nav, resequence_nav, collections_nav, collection_name,
            collection_selector_list, run_list]

@app.callback(
    [Output("run-list", "options"), 
     Output("collection-selector", "options"),
     Output("group-list", "options"),
     Output("species-list", "options")],
    [Input("form-species-source", "value"),
    Input("selected-collection", "data")
    ]
)
@cache.memoize(timeout=cache_timeout)  # in seconds
def update_run_options(form_species, selected_collection):
    return filter_update_run_options(form_species, selected_collection)


@app.callback(
    Output("collection-selector", "value"),
    [Input("selected-collection", "data")]
)
def update_selected_collection(selected_collection):
    return selected_collection

@app.callback(
    [Output("run-list", "value"),
     Output("group-list", "value"),
     Output("species-list", "value"),
     Output("qc-list", "value"),
     Output("samples-form", "value")],
    [Input("param-store", "data")]
)
def update_filter_values(param_store):
    return filter_update_filter_values(param_store)


@app.callback(
    Output("collection-link", "href"),
    [Input("collection-selector", "value")],
    [State("url", "pathname")]
)
def update_collection_button_f(collection, pathname):
    return update_collection_button(collection, pathname)

@app.callback(
    Output("page-n",
            "children"),
    [Input("prevpage", "n_clicks_timestamp"),
        Input("prevpage2", "n_clicks_timestamp"),
        Input("nextpage", "n_clicks_timestamp"),
        Input("nextpage2", "n_clicks_timestamp")],
    [State("page-n", "children"),
        State("max-page", "children")]
)
def next_page(prev_ts, prev_ts2, next_ts, next_ts2, page_n, max_page):
    return samples_next_page(prev_ts, prev_ts2, next_ts, next_ts2, page_n, max_page)


@app.callback(
    Output("sample-report", "children"),
    [Input("page-n", "children"),
    Input("sample-store", "data")]
        )
def fill_sample_report(page_n, sample_store):
    page_n = int(page_n)
    sample_ids = list(
        map(lambda x: x["_id"], sample_store))
    if len(sample_ids) == 0:
        return None
    
    data_table = import_data.filter_all(
        sample_ids=sample_ids,
        pagination={"page_size": SAMPLE_PAGESIZE, "current_page": page_n})
    max_page = len(sample_store) // SAMPLE_PAGESIZE
    # We need to have fake radio buttons with the same ids to account for times 
    # when not all SAMPLE_PAGESIZE samples are shown and are not taking the ids required by the callback
    html_fake_radio_buttons = html.Div([dcc.RadioItems(
        options=[
            {'label': '', 'value': 'nosample'}
        ],
        value='noaction',
        id="sample-radio-{}".format(n_sample)
    ) for n_sample in range(len(data_table), SAMPLE_PAGESIZE)], style={"display": "none"})
    return [
        html.H4("Page {} of {}".format(page_n + 1, max_page + 1)),
        html.Div(children_sample_list_report(data_table)),
        html_fake_radio_buttons,
        admin.html_qc_expert_form(),
        html.H4("Page {} of {}".format(page_n + 1, max_page + 1)),
        dcc.ConfirmDialog(
            id='qc-confirm',
            message='Are you sure you want to send sample feedback?',
        )
    ]


@app.callback(
    Output("sample-store", "data"),
    [Input("apply-filter-button", "n_clicks"),
     Input("param-store", "data"),
     Input("selected-collection", "data")],
    [State("run-list", "value"),
     State("species-list", "value"),
     State("form-species-source", "value"),
     State("group-list", "value"),
     State("qc-list", "value"),
     State("samples-form", "value"),
     State("sample-store", "data"),
     State("date-sequenced", "start_date"),
     State("date-sequenced", "end_date"),
     ]
)
def update_selected_samples(n_clicks, param_store, collection_name,
                            run_names, species_list,
                            species_source, group_list, qc_list,
                            sample_names, prev_sample_store,
                            date_seq_start, date_seq_end):
    date_range = [date_seq_start, date_seq_end]
    for i in range(2):
        if date_range[i] is not None:
            date_range[i] = datetime.datetime.strptime(re.split('T| ', date_range[i])[0], '%Y-%m-%d')

    if sample_names is not None and sample_names != "":
        sample_names = sample_names.split("\n")
    else:
        sample_names = param_store.get("sample_names", [])
    if not run_names:
        run_names = param_store.get("run", [])
    if not group_list:
        group_list = param_store.get("group", [])
    if not species_list:
        species_list = param_store.get("species", [])
    if not qc_list:
        qc_list = param_store.get("qc", [])
    if not date_range[0]:
        date_range[0] = param_store.get("date_seq_start", None)
    if not date_range[1]:
        date_range[1] = param_store.get("date_seq_end", None)

    #override if selected collection
    if collection_name is not None:
        run_names = [collection_name]

    if (date_range[0] is None and
            date_range[1] is None):
        date_range = None

    if (n_clicks == 0 and
        sample_names == [] and
        run_names == [] and
        group_list == [] and
        species_list == [] and
        qc_list == [] and
        date_range is None):
        samples = prev_sample_store
    else:
        

        samples = import_data.filter_all(
            species=species_list, species_source=species_source,
            group=group_list, qc_list=qc_list,
            run_names=run_names,
            sample_names=sample_names,
            date_range=date_range,
            projection={"name": 1})

        if "_id" in samples:
            samples["_id"] = samples["_id"].astype(str)
        samples = samples.to_dict('records')
    # if deleted_samples:
    #     samples = [s for s in samples if s["_id"] not in deleted_samples]
    return samples


@app.callback(
    [
        Output("filter-sample-count", "children"),
        Output("datatable-ssi_stamper", "data"),
        Output("datatable-ssi_stamper", "virtualization")
    ],
    [
     Input("placeholder0", "children"),
     Input("sample-store", "data")
    ],
)
def update_filter_table(_, sample_store):
    if len(sample_store) == 0:
        return ["0", [{}], False]
    sample_ids = list(
        map(lambda x: x["_id"], sample_store))

    samples = import_data.filter_all(
        sample_ids=sample_ids)

    samples = generate_table(samples)
    if len(sample_store) > 500:
        virtualization = True
    else:
        virtualization = False
    return [len(sample_store), samples.to_dict("rows"), virtualization]

@app.callback(
    Output("tsv-download", "children"),
    [Input("generate-download-button", "n_clicks")],
    [State("run-list", "value"),
     State("species-list", "value"),
     State("form-species-source", "value"),
     State("group-list", "value"),
     State("qc-list", "value"),
     State("samples-form", "value")],
    prevent_initial_call=True
)
def generate_download_button(download_button,
                            run_names, species_list,
                            species_source, group_list, qc_list,
                            sample_names):
    if download_button == 0:
        return None
    else:
        if sample_names is not None and sample_names != "":
            sample_names = sample_names.split("\n")

        tests_df = import_data.filter_all(species=species_list, species_source=species_source,
                                                       group=group_list, qc_list=qc_list,
                                                       run_names=run_names,
                                                       sample_names=sample_names,
                                                       pagination=None)
    # return samples.to_dict()
    if not len(tests_df):
        return None

    tests_df = generate_table(tests_df)

    rename_dict = {item["id"]: item["name"]
        for item in global_vars.COLUMNS}

    renamed = tests_df.rename(rename_dict, axis='columns')

    missing_columns = [a for a in list(
        rename_dict.values()) if not a in list(renamed.columns)]

    # add missing columns
    for column in missing_columns:
        renamed[column] = np.nan

    # reorder columns
    renamed = renamed[list(rename_dict.values())]

    csv_string_eur = renamed.to_csv(
        index=False, encoding="utf-8", sep=";", decimal=",")
    tsv_string_us = renamed.to_csv(index=False, encoding="utf-8", sep="\t")
    full_csv_string_eur = 'data:text/csv;charset=utf-8,' + \
        urlparse.quote(csv_string_eur)
    full_tsv_string_us = 'data:text/tab-separated-values;charset=utf-8,' + \
        urlparse.quote(tsv_string_us)
    return [
        html.A("(tsv, US format)",
                href=full_tsv_string_us,
                download='report.tsv'),
        " - ",
        html.A("(csv, EUR Excel format)",
                href=full_csv_string_eur,
                download='report.csv')
    ]


@app.callback(
    [Output("plot-species", "value"),
     Output("plot-species", "options")],
    [Input("sample-store", "data"),
     Input("plot-species-source", "value")],
    [State("plot-species", "value")]
)
def aggregate_species_dropdown_f(sample_store, plot_species, selected_species):
    return aggregate_species_dropdown(sample_store, plot_species, selected_species)


@app.callback(
    [Output("pipeline-table", "data"),
     Output("pipeline-table", "columns"),
     Output("pipeline-table", "style_data_conditional"),
     Output("rerun-samples", "options"),
     Output("rerun-components", "options")],
    [Input("sample-store", "data"),
     Input("table-interval", "n_intervals")]
)
def pipeline_report_data_f(sample_store, _):
    return pipeline_report_data(sample_store)


@app.callback(
    [Output("summary-plot", "figure"),
     Output("mlst-plot", "figure")],
    [Input("plot-species", "value")],
    [State("sample-store", "data"),
    State("plot-species-source", "value")]
)
@cache.memoize(timeout=cache_timeout)  # in seconds
def update_aggregate_fig_f(selected_species, samples, plot_species_source):
    return update_aggregate_fig(selected_species, samples, plot_species_source)

@app.callback(Output("pipeline-rerun", "data"),
              [Input("pipeline-table", "active_cell"),
               Input("pipeline-table", "derived_viewport_data"),
               Input("rerun-add-components", "n_clicks"),
               Input("rerun-add-samples", "n_clicks"),
               Input("rerun-add-failed", "n_clicks")], 
              [State("pipeline-table", "columns"),
               State("pipeline-rerun", "derived_viewport_data"),
               State("rerun-components", "value"),
               State("rerun-samples", "value")])
def update_rerun_table_f(active, table_data, n_click_comp, n_click_samp,
                       n_click_fail, columns, prev_data, rerun_comp,
                       rerun_samp):
    return update_rerun_table(active, table_data, n_click_comp, n_click_samp,
                              n_click_fail, columns, prev_data, rerun_comp,
                              rerun_samp)


@app.callback(
    [Output("rerun-output", "children"),
     Output("rerun-output", "is_open")],
    [Input("rerun-button", "n_clicks")],
    [State("pipeline-rerun", "derived_viewport_data")],
    prevent_initial_call=True
)
def rerun_components_button_f(n_clicks, data):
    return rerun_components_button(n_clicks, data, config["rerun"])


@app.callback(Output('qc-confirm', 'displayed'),
              [Input('feedback-button', 'n_clicks_timestamp')],
              prevent_initial_call=True)
def display_confirm_feedback(button):
    if button is not None:
        return True
    return False


@app.callback(
    Output("qc-feedback", "children"),
    [Input("qc-confirm", "submit_n_clicks")],
    [State("qc-user-1", "value")] + [State("sample-radio-{}".format(n), "value")
                                     for n in range(SAMPLE_PAGESIZE)] +
                                    [State("sample_reason-{}".format(n), "value")
                                     for n in range(SAMPLE_PAGESIZE)],
    prevent_initial_call=True
)
def submit_user_feedback(_, user, *args):
    if (config["feedback_enabled"]):
        feedback_pairs = []
        for i in range(int(len(args)/2)):
            val = args[i]
            reason = args[int(len(args)/2) + i]
            if val != "noaction":
                if val.startswith("OK_"):
                    feedback_pairs.append((val[3:], "OK", reason))
                elif val.startswith("CF_"):
                    feedback_pairs.append((val[3:], "resequence", reason))
                elif val.startswith("OT_"):
                    feedback_pairs.append((val[3:], "other", reason))
        if len(feedback_pairs) > 0:
            email_config = {
                "email_from": config["email_from"],
                "email_to": config["email_to"]
            }
            import_data.add_batch_user_feedback_and_mail(feedback_pairs, user, email_config)
            return "Feedback saved"
    return []


@app.callback(
    Output("link-to-files-div", "children"),
    [Input("sample-store", "data")],
)
def link_to_files_f(data):
    return link_to_files(data)


app.clientside_callback(
    ClientsideFunction(
        namespace='clientside',
        function_name='enable_selection_button'
    ),
    [Output("add-collection-button", "disabled"),
     Output("remove-collection-button", "disabled")],
    [Input("datatable-ssi_stamper", "derived_viewport_selected_rows")]
)

@app.callback(
    Output("selection-modal-status", "children"),
    [Input("selection-modal-create", "n_clicks")],
    [State("datatable-ssi_stamper", "derived_viewport_selected_rows"),
     State("datatable-ssi_stamper", "derived_viewport_data"),
     State("virtual-run-selector", "value"),
     State("virtual_run_name", "value")
    ],
    prevent_initial_call=True
)
def add_selected_to_collection(n_clicks, selected_rows, row_data, selected_virtual, name):
    data = []
    message = ""
    if (selected_rows is not None and len(selected_rows)):
        for e in selected_rows:
            row = {
                "name": row_data[e]["name"],
                "_id": row_data[e]["id"]
            }
            if (row not in data):
                data.append(row);
    if n_clicks > 0:
        if selected_virtual == "_new":
            message = virtual_runs.create_run(data, name)
        else:
            message = virtual_runs.add_samples_to_virtual_run(data, selected_virtual)

    return message

@app.callback(
    Output("remove-modal-status", "children"),
    [Input("remove-collection-modal-remove", "n_clicks")],
    [State("datatable-ssi_stamper", "derived_viewport_selected_rows"),
     State("datatable-ssi_stamper", "derived_viewport_data"),
     State("selected-collection", "data")
    ],
    prevent_initial_call=True
)
def remove_selected_from_collection(n_clicks, selected_rows, row_data, name):
    data = []
    message = ""
    if (selected_rows is not None and len(selected_rows)):
        for e in selected_rows:
            if (row_data[e]["id"] not in data):
                data.append(row_data[e]["id"]);
    if n_clicks > 0:
        message = virtual_runs.remove_samples_from_virtual_run(data, name)

    return message


# app.clientside_callback(
#     ClientsideFunction(
#         namespace='clientside',
#         function_name='open_close_selection_modal'
#     ),
#     [Output("selection-modal", "is_open"), Output("selection-modal-body", "children")],
#     [Input("selection-modal-open", "n_clicks"), Input("selection-modal-close", "n_clicks")],
#     [State("selection-modal", "is_open"),
#      State("selected-samples", "data"),
#     ]
# )
@app.callback(
    Output("add-collection-modal", "is_open"),
    [Input("add-collection-button", "n_clicks"), Input("add-collection-modal-close", "n_clicks")],
    [State("add-collection-modal", "is_open")]
)
def open_close_add_modal(n1, n2, is_open):
    trigger = dash.callback_context.triggered[0]

    if (n1 or n2) and trigger["value"] is not None:
        return not is_open
    return is_open

@app.callback(
    Output("remove-collection-modal", "is_open"),
    [Input("remove-collection-button", "n_clicks"), Input("remove-collection-modal-close", "n_clicks")],
    [State("remove-collection-modal", "is_open")]
)
def open_close_remove_modal(n1, n2, is_open):
    if (n1 or n2):
        return not is_open
    return is_open

app.clientside_callback(
    ClientsideFunction(
        namespace='clientside',
        function_name='fill_full_virtual_name'
    ),
    Output("virtual_run_full_name", "children"),
    [Input("virtual_run_name", "value")],
)

server = app.server # Required for gunicorn

def main_debug():
    app.run_server(debug=True, host="0.0.0.0", dev_tools_hot_reload=True)

if __name__ == '__main__':
    # 0.0.0.0 exposes the app to the network.
    main_debug()
