import pandas as pd
import dash_html_components as html

def html_table(data, **kwargs):
    rows = []
    for data_row in data:
        if isinstance(data_row, list):
            rows.append(
                html.Tr([html.Td(data_cell, className="cell") for data_cell in data_row], className="trow"))
        else:
            rows.append(html.Tr([html.Td(data_cell, className="cell") for data_cell in data_row["list"]],
                               className=data_row["className"] + " trow"))
    return html.Table(rows, className="bifrost-table", **kwargs)


def html_td_percentage(value, color):
    string = str(round(float(value) * 100, 2)) + "%"
    if pd.isna(value):
        fill = "0%"
    else:
        fill = string
    return html.Td(
        html.Div([
            html.Span(string, className="val"),
            html.Span(
                className="bar",
                style={"backgroundColor": color, "width": fill}
            )
        ], className="wrapper"
        ), className="data-colored cell")
