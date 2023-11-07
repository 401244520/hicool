from dash import Dash, html, dcc, Input, Output, callback
import pandas as pd
import plotly.express as px
import glob
from scipy.interpolate import interp1d
import numpy as np

def track_handler(track_path,chrom,start,end,bin_num=1000):
    suffix = track_path.split(".")[-1]
    if suffix in ["bw","BigWig","bigwig"]:
        import pyBigWig
        with pyBigWig.open(track_path) as bw_file:
            coverage = bw_file.stats(chrom,start,end,type="mean",nBins=bin_num)
    elif suffix in ["bdg","BedGraph","bedgraph"]:   
        df = pd.read_csv(track_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'value'])
        df_region = df[(df['chrom'] == chrom) & (df['end'] >= start) & (df['start'] <= end)]
        coverage = np.zeros(end - start)
        for i,interval in df_region.iterrows():
            coverage[interval.start - start:interval.end - start] = interval.value
    else :
        raise ValueError(f"Do not support {suffix} format!")
    if bin_num is None:
        bin_num = (end - start) // 1000 # 默认为1kb分辨率
    return resample_vec(coverage,bin_num)

def bw_track():
    bw_list = glob.glob("../../../tests/data/*[bw,bdg]")
    return html.Div([
        html.Div([
            dcc.Dropdown(
                    bw_list,
                    bw_list[0],
                    id='bw_path',
                    ),
            html.Div([
                dcc.Dropdown(id='chrom', style={'height': '36px'})
            ], style={'display': 'inline-block', 'width': '10%', 'height': '22px'}),
            html.Div([
                dcc.Input(id='start', placeholder='start:0',
                          type='int', min=0, value=0, style={'height': '30px'}),
                dcc.Input(id='end', placeholder='end:',
                          type='int', style={'height': '30px'}),
                html.Button(id='submit_button',
                            children='Show Heatmap', style={'height': '30px'}),
            ], style={'display': 'inline-block'}),
        ]),
        html.Div([
            dcc.Slider(
                id='bin_num', min=100, max=1000000,
                ),
            dcc.Graph(
                id='bw_graph',
            )
        ]),
    ])


@callback(
    Output('bw_graph', 'figure'),
    Input('bw_path', 'value'),
    Input('chrom', 'value'),
    Input('start', 'value'),
    Input('end', 'value'),
    Input('bin_num', 'value')
    )
def update_graph(bw_path,chrom,start,end,bin_num):
    x = np.arange(start,end,bin_num)
    y = track_handler(bw_path,chrom,start,end,bin_num)
    fig = px.line(x=x,y=y,title=f"{chrom}:{start}-{end}")
    fig.update_layout(
        margin={'l': 40, 'b': 40, 't': 10, 'r': 0}, hovermode='closest')
    return fig
