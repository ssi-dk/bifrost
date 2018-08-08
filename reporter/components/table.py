
import dash_html_components as html

def html_table(data, **kwargs):
    rows = []
    for data_row in data:
        if isinstance(data_row, list):
            rows.append(html.Tr([html.Td(data_cell) for data_cell in data_row]))
        else:
            print(data_row)
            rows.append(html.Tr([html.Td(data_cell) for data_cell in data_row["list"]],
                               className=data_row["className"]))
    return html.Table(rows, **kwargs)


def html_td_percentage(value, color):
    string = str(round(float(value) * 100, 2)) + "%"
    return html.Td(
        html.Div([
            html.Span(string, className="val"),
            html.Span(
                className="bar",
                style={"backgroundColor": color, "width": string}
            )
        ], className="wrapper"
        ), className="data-colored")
