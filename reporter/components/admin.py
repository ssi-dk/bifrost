import os
import dash_html_components as html
import dash_core_components as dcc


def html_qc_expert_form():
    if ("REPORTER_ADMIN" in os.environ and os.environ["REPORTER_ADMIN"] == "True"):
        return html.Div([
            html.Div(html.H5("Submit Supplying Lab Feedback:",
                             className="nomargin"), className="six columns"),
            html.Div(dcc.Input(id="qc-user-1", placeholder="Your initials (e.g. MBAS)",
                               value="", type="text"), className="three columns"),
            html.Div(html.Button("Submit QC values",
                                 id="feedback-button"), className="three columns")
        ], className="row mb-1")
    else:
        return None


def sample_radio_feedback(sample, n_sample):
    if ("REPORTER_ADMIN" in os.environ and os.environ["REPORTER_ADMIN"] == "True"):
        return (html.Div([
            html.Div(html.Label("Supplying Lab Feedback:"),
                     className="three columns"),
            html.Div(
                dcc.RadioItems(
                    options=[
                        {'label': 'Accept', 'value': 'A_{}'.format(
                                            sample["_id"])},
                        {'label': 'Resequence', 'value': 'R_{}'.format(
                            sample["_id"])},
                        # {'label': 'Other', 'value': 'O_{}'.format(
                        #     sample["_id"])},
                        {'label': 'No action', 'value': 'noaction'}
                    ],
                    value='noaction',
                    id="sample-radio-{}".format(n_sample),
                    labelStyle={'display': 'inline-block'}), className="nine columns")
        ], className="row mt-1"))
    else:
        return None
