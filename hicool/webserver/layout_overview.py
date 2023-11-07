import glob
import plotly.express as px
import pandas as pd
from dash import dash_table, html, dcc, Input, Output, callback, State

try:
    from hicool.api import HiCool
except:
    import sys
    sys.path.append("/home/wzl/Workspace/HiCool")
    from hicool.api import HiCool


def hicool_view(datasets):
    hicool_list = glob.glob(f"{datasets}/uploads/*hicool")
    print(str(hicool_list) + " \n \n")
    return html.Div([
        html.Div([
            dcc.Dropdown(
                options=[{'label': file, 'value': file}
                    for file in hicool_list],
                value=hicool_list[0],
                id='hicool_path',
            ),
            dcc.RadioItems(id='hicool_content', inline=True),
            dcc.RadioItems(id='hicool_bottom', inline=True),
        ]),
        html.Div(id='hicool_info'),
    ])


@callback(
    Output('hicool_content', 'options'),
    Output('hicool_content', 'value'),
    Input('hicool_path', 'value'))
def update_content(hicool_path):
    global hc
    hc = HiCool(hicool_path)
    content = sorted(list(hc.show()))
    return content, content[0]


@callback(
    Output('hicool_bottom', 'options'),
    Output('hicool_bottom', 'value'),
    Input('hicool_path', 'value'),
    Input('hicool_content', 'value'))
def update_botton(hicool_path, hicool_content):
    # hc = HiCool(hicool_path)
    bottom = sorted(list(hc.show(hicool_content)))
    return bottom, bottom[0]


@callback(
    Output('hicool_info', 'children'),
    Input('hicool_path', 'value'),
    Input('hicool_content', 'value'),
    Input('hicool_bottom', 'value'))
def update_view(hicool_path, hicool_content, hicool_bottom):
    # hc = HiCool(hicool_path)
    root = hicool_content+"/"+hicool_bottom
    data = hc.show(root)
    if type(data) == set:
        return html.Div([
            html.H6(
                f'{root} is hdf5.group, including {len(data)} groups or datasets.'
            ),
            html.H6("You can search group information in the window below."),
            html.Div([
                "hdf5 root: ",
                dcc.Input(id='root_input', value=root +
                          "/", type='text', size=30),
                html.Button(id='submit_button', children='Show Contents'),
            ]),
            html.Br(),
            html.Div(id='root_output'),
        ])
    else:
        return [show_data(data, root)]


@callback(
    Output('root_output', 'children'),
    Input('hicool_path', 'value'),
    State('root_input', 'value'),
    Input('submit_button', 'n_clicks'),
)
def update_output_div(hicool_path, root_input, n_clicks):
    # hc = HiCool(hicool_path)
    try:
        data = hc.show(root_input)
        if type(data) == set:
            index = sorted(list(data))
            if len(index) > 10 :
                index = index[:10] + ["......"]
            return f'{root_input} : {index}'
        else:
            return [show_data(data, root_input)]
    except:
        return 'Please enter valid path.'


def show_data(data, root_path):
    root = root_path.split("/")[-1]
    # print(df[root].dtypes)
    try:
        df = pd.DataFrame(data, columns=[root])
        return html.Div([
            html.Div([
                html.H6(
                    f"{root_path} including {len(df)} values in {df[root].dtypes} format."),
                html.H6(
                    f"Ranging {df[root].round(2).min()} to {df[root].round(2).max()} with {df[root].round(2).median()} median."),
                dash_table.DataTable(
                    id='table',
                    columns=[{"name": i, "id": i}
                             for i in df.columns],
                    data=df.to_dict('records'),
                    page_size=10
                )
            ], style={'width': '30%', 'display': 'inline-block'}),
            html.Div([
                dcc.Graph(figure=px.histogram(df,nbins=100,histnorm='percent'))
            ], style={'width': '60%', 'display': 'inline-block',"height":"400px"}),
        ])
    except:
        df = pd.DataFrame(data)
        for col in df.columns:
            if df[col].dtypes in [object, '|S5']:
                df[col] = df[col].apply(lambda x: x.decode())
        df.columns = df.columns.astype(
            str) + ":" + df.dtypes.apply(lambda x: str(x))
        return html.Div([
            html.H6(
                f"{root_path} shapes {df.shape} , col name shows dtype (max show 10 cols) ."),
            dash_table.DataTable(
                id='table',
                columns=[{"name": i, "id": i}
                         for i in df.columns[:10]],
                data=df.iloc[:, :10].to_dict('records'),
                page_size=10
            )
        ], style={'width': '30%', 'display': 'inline-block'})
