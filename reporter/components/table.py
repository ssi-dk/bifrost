
import dash_html_components as html

def html_table(data, **kwargs):
    return html.Table([
        html.Tr([html.Td(data_cell) for data_cell in data_row])
        for data_row in data
    ], **kwargs)


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
