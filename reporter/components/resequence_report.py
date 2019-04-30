
import dash_html_components as html


import components.mongo_interface
import components.import_data as import_data


def resequence_report(run_name):
    run_checker_link = html.H4(html.A(
        "Run Checker Report", href="/{}".format(run)))

    last_runs = import_data.get_last_runs(run["name"], 12) #Get last 12 runs
    last_runs_names = [run["name"] for run in last_runs]
    prev_runs_dict = import_data.get_sample_QC_status(last_runs)
    header = html.Tr([html.Th(html.Div(html.Strong("Sample")),
                                className="rotate")] +
                        list(map(lambda x: html.Th(html.Div(html.Strong(x)),
                                                className="rotate"),
                                last_runs_names)))
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
                    # to account for libray fails
                    sample_all_OKs = False
                row.append(
                    html.Td(status,
                            className="center status-" + className))
            else:
                row.append(html.Td("-", className="center status-0"))

        if not sample_all_OKs:
            rows.append(html.Tr(row))
    table = html.Table(rows, className="unset-width-table")
    update_notice += (" SL: Supplying Lab, CF: Core Facility, CF(LF): "
                        "Core Facility (Library Fail). -: No data. "
                        "*: user submitted")
    return [
        run_checker_link,
        html.P(update_notice),
        table
    ]
