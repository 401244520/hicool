from dash import Dash, html, dcc, Input, Output, callback,State
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import glob
import numpy as np
import pyBigWig


def bw_track(datasets):
    bw_list = glob.glob(f"{datasets}/uploads/*bw")
    return html.Div([
        html.Div([
            html.Div([
                dcc.Dropdown(id='bw_path',
                    options=bw_list,
                    value=bw_list[0],
                    style={'height': '28px'},
                )
            ],style={'display': 'inline-block','textAlign': 'center', 'width': '80%'}),
            html.Div([
                html.Button(id='add_button',n_clicks=0,
                            children='Add BigWig', style={'height': '30px'}),
            ], style={'display': 'inline-block', 'width': '10%','position': 'relative', 'bottom': '5px'}
            )
        ]),
        html.Div([
            html.Div([
                dcc.Dropdown(id='chrom', style={'height': '28px'}),
            ],style={'display': 'inline-block', 'width': '10%','textAlign': 'center'}),
            html.Div([
                dcc.Input(id='start', placeholder='start:0',
                          type='number', min=0, value=0, style={'height': '30px','textAlign': 'center'}),
                dcc.Input(id='end', placeholder='end:',
                          type='number',min=0, style={'height': '30px','textAlign': 'center'}),
                html.Button(id='submit_button',n_clicks=0,
                            children='Show BigWig', style={'height': '30px'}),
            ], style={'display': 'inline-block','position': 'relative', 'bottom': '5px'}),
        ]),
        html.Div([
            dcc.RadioItems(id='bw_type', value='mean', inline=True,
                            options=['mean','max','min','coverage']),
            dcc.Checklist(id='bw_list',options=[],value=[])
        ]),
        html.Div([
            dcc.Slider(id='bin_num', min=10, max=10000,value=1000),
        html.Div([
            dcc.Graph(id='bw_graph'),
        ], style={'display': 'inline-block', 'width': '90%'}),
        html.Div([
            dcc.Slider(id='bw_max', min=0, max=10, vertical=True),
        ], style={'height': '-50px', 'display': 'inline-block', 'width': '5%'}),
        html.Div([
            dcc.Slider(id='bw_min', min=0, max=10, vertical=True),
        ], style={'height': '-50px', 'display': 'inline-block', 'width': '5%'})
        ]),
    ])


@callback(
    Output('chrom', 'options'),
    Output('chrom', 'value'),
    Output('end', 'value'),
    Input('bw_path', 'value'))
def select_bw(bw_path):
    with pyBigWig.open(bw_path) as bw_file:
        chroms = list(bw_file.chroms().keys())
        end = bw_file.chroms(chroms[0])
    return chroms,chroms[0],end

@callback(
    Output('bw_list', 'options'),
    State('bw_list', 'options'),
    State('bw_path', 'value'),
    Input('add_button', 'n_clicks'))
def add_bw(bw_list,bw_path, n_clicks):
    if (bw_path not in bw_list) and n_clicks > 0:
        bw_list.append(bw_path)
    return bw_list


@callback(
    Output('bw_graph', 'figure'),
    Input('bw_list', 'value'),
    State('bw_path', 'value'),
    Input('chrom', 'value'),
    State('start', 'value'),
    State('end', 'value'),
    Input('bin_num', 'value'),
    Input('bw_type', 'value'),
    Input('bw_min', 'value'),
    Input('bw_max', 'value'),
    Input('submit_button', 'n_clicks'),
    )
def update_graph(bw_list,bw_path,chrom,start:int,end:int,bin_num,bw_type,bw_min,bw_max,n_clicks):
    print(chrom,start,end,bin_num,bw_list)
    if bw_list == []:
        bw_list = [bw_path]
    fig = make_subplots(rows=len(bw_list), cols=1) 
    for i,bw_path in enumerate(bw_list):
        with pyBigWig.open(bw_path) as bw_file:
            chrom_len = bw_file.chroms(chrom)
            if (end is None) or (end > chrom_len):
                end = chrom_len
            x = np.linspace(start,end,bin_num)
            y = bw_file.stats(chrom,start,end,type=bw_type,nBins=bin_num)
        name = bw_path.split("/")[-1].split(".")[0]
        trace = go.Scatter(x=x, y=y, mode='lines', fill='tozeroy',name=name)
        fig.add_trace(trace, row=i+1, col=1)
        fig.update_yaxes(title_text = bw_type, range=[bw_min,bw_max], row=i+1, col=1)
        fig.update_xaxes(title_text = name, row=i+1, col=1)
    fig.update_layout(
        title=f"{chrom}:{start}-{end}", title_x=0.5,
        height = len(bw_list)*200 + 200,
        margin={'l': 40, 'b': 100, 't': 100, 'r': 40})
    return fig
