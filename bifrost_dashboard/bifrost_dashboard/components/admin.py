import os
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc


def html_qc_expert_form():
    if ("REPORTER_ADMIN" in os.environ and os.environ["REPORTER_ADMIN"] == "True"):
        return html.Div([
            dbc.Row([
                dbc.Col(html.H5("Submit Supplying Lab Feedback:",
                                className="nomargin"), width=6),
                dbc.Col(dcc.Input(id="qc-user-1", placeholder="Your initials (e.g. MBAS)",
                                  value="", type="text"), width=3),
                html.Div(id="placeholder0", style={"display": "none"}),
                dbc.Col([
                    html.Button("Submit QC values",
                                id="feedback-button"),
                    html.Div(id="qc-feedback", className="qc-feedback")], width=3)
            ]),
            html.Div([
                html.Details([
                    html.Summary("How it works (click to show)"),
                    html.P("""If you want to send a sample to resequence,
                    go to the sample card and select "Resequence" in 
                    the Suplying Lab Feedback section. The admin will get
                    an email with the change. You can also "Accept" samples
                    with the "Supplying lab" warning and their data will
                    be used in the future to adjust QC thresholds.""")
                ])
            ], className="mb-1")

        ])
        
    else:
        return None


def sample_radio_feedback(sample, n_sample):
    if ("REPORTER_ADMIN" in os.environ and os.environ["REPORTER_ADMIN"] == "True"):
        return (html.Div([
            html.Label("Supplying Lab Feedback:", className="mr-4",
                       style={'display': 'inline-block'}),
            dcc.RadioItems(
                options=[
                    {'label': ' Accept', 'value': 'OK_{}'.format(
                                        sample["_id"])},
                    {'label': ' Resequence', 'value': 'CF_{}'.format(
                        sample["_id"])},
                    {'label': ' Other', 'value': 'O_{}'.format(
                        sample["_id"])},
                    {'label': ' No action', 'value': 'noaction'}
                ],
                value='noaction',
                id="sample-radio-{}".format(n_sample),
                labelStyle={
                    'margin': '0 0.5rem 0.5rem 0'},
                style={'display': 'inline-block'}),
        ]))
    else:
        return None
