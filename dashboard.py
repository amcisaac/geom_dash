# Import packages
from dash import Dash, html, dash_table, dcc, callback, Output, Input
import pandas as pd
import plotly.express as px
import json
import plotly.graph_objects as go

with open('angles_test2_rnmm.json','r') as jsonfile:
    test = dict(json.load(jsonfile))

test_keys = list(test.keys())
test_df = pd.DataFrame(test[test_keys[0]])

@callback(
    Output(component_id = 'controls-and-graph',component_property='figure'),
    Input(component_id = 'param-id',component_property='value')
)
def make_figure(key):
    df = pd.DataFrame(test[key])
    layout = go.Layout(
        title=df['ident'][0] +': '+ key,
        xaxis=dict(
            title="MM value"
        ),
        yaxis=dict(
            title="QM value"
    ) )
    fig = go.Figure(layout=layout)
    fig.add_trace(go.Scatter(x = [df['qm_values'].min(),df['qm_values'].max()],y = [df['qm_values'].min(),df['qm_values'].max()],mode='lines',line=dict(color='black',dash='dash'),name='x=y'))

    fig.add_trace(go.Scatter(x = df['espaloma_values'], y = df['qm_values'],mode='markers',marker=dict(color='blue')))
    fig.add_trace(go.Scatter(x = [df['sage_value'][0]], y = [df['sage_value'][0]],mode='markers', marker=dict(size=[15],color='red'),marker_symbol='star',name='Sage'))

    return fig


# @callback(
#     [
#         Output("graph-container", "children", allow_duplicate=True),
#         Output("smirks_input", "value", allow_duplicate=True),
#     ],
#     Input("previous", "n_clicks"),
#     prevent_initial_call=True,
# )
# def previous_button(_):
#     global CUR_SMIRK
#     if CUR_SMIRK >= 1:
#         CUR_SMIRK -= 1
#     return make_fig(cur_record(), NCLUSTERS), SMIRKS[CUR_SMIRK]
#
#
# @callback(
#     [
#         Output("graph-container", "children", allow_duplicate=True),
#         Output("smirks_input", "value", allow_duplicate=True),
#     ],
#     Input("next", "n_clicks"),
#     prevent_initial_call=True,
# )
# def next_button(_):
#     global CUR_SMIRK
#     if CUR_SMIRK < len(SMIRKS) - 1:
#         CUR_SMIRK += 1
#     return make_fig(cur_record(), NCLUSTERS), SMIRKS[CUR_SMIRK]





app = Dash(__name__)

app.layout = html.Div([
    html.Div(children='QM vs MM geometry parameters'),

    dcc.Graph(figure=make_figure(test_keys[0]), id = 'controls-and-graph'),
    html.Div([

        html.Div([
            dcc.Dropdown(
                test_keys,
                'Parameter',
                id='param-id'
            )
        ])
    ])
    # html.Button("Previous", id="previous", n_clicks=0),
    # html.Button("Next", id="next", n_clicks=0),

])

if __name__ == '__main__':
    app.run(debug=True)
