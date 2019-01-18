import os
import dash_html_components as html

def selected_samples_div():
    if ADMIN == True:
        return (html.Div([
            "Selected samples: ",
            html.Button(
                "OK", className="button passfail", id="qc-pass-button"),
            html.Button("suppl. lab", className="button passfail",
                        id="qc-sl-button"),
            html.Button("core fac.", className="button passfail",
                        id="qc-cf-button")
        ], className="u-pull-right", id="qc-buttons"))
    else:
        return None
