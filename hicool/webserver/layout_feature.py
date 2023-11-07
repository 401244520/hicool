from dash import Dash, html, dcc, Input, Output, callback,State
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import glob
import numpy as np

try:
    from hicool.tools import HiCool
    from hicool.function.dataloader import load_cool_region
except:
    import sys
    sys.path.append("/home/wzl/Workspace/HiCool")
    from hicool.tools import HiCool
    from hicool.function.dataloader import load_cool_region

def hicool_features(datasets):
    hicool_list = glob.glob(f"{datasets}/uploads/*hicool")
    return html.Div([
        html.Div([
            dcc.Dropdown(id='hicool_path',
                    options=hicool_list,
                    value=hicool_list[0],
            )
        ],style={'display': 'inline-block','textAlign': 'center', 'width': '50%'}),
        
        html.Div([
            dcc.Dropdown(id='features', multi=True)
        ], style={'display': 'inline-block','textAlign': 'center', 'width': '50%'}),
        
        html.Div([
            html.Div([
                dcc.Dropdown(id='feature_chrom')
            ], style={'display': 'inline-block','textAlign': 'center', 'height': '22px', 'width': '12%'}),
            html.Div([
                dcc.Input(id='feature_start', placeholder='start:0', 
                        type='int', min=0, value=0, style={'height': '30px','textAlign': 'center'}),        
            ], style={'display': 'inline-block', 'width': '10%','margin': '0 1%'}),
            html.Div([
                dcc.Input(id='feature_end', placeholder='end:', type='int',
                          style={'height': '30px','textAlign': 'center'}),
            ], style={'display': 'inline-block', 'width': '10%', 'margin': '0 1%'}),
            html.Div([
                html.Button(id='feature_button', children='Show Features', 
                            style={'height': '30px','textAlign': 'center'}),
            ], style={'display': 'inline-block', 'width': '10%', 'margin': '0 1%'}),
        ]),
                
        html.Div([
            dcc.RadioItems(id='hicool_field1', value='count')
        ]),
        html.Div([
           dcc.Slider(id='hicool_zmax1', min=0, max=100),
        ], style={'display': 'inline-block', 'width': '50%'}),
        html.Div([
            dcc.Graph(id='compare_heatmap'),
        ], style={'display': 'inline-block', 'width': '80%'}),
    ])

@callback(
    Output('hicool_cell1', 'options'),
    Output('hicool_cell1', 'value'),
    Output('hicool_cell2', 'options'),
    Output('hicool_cell2', 'value'),
    Output('hicool_field1', 'options'),
    Output('hicool_features', 'options'),
    Output('hicool_chrom1', 'options'),
    Output('hicool_chrom1', 'value'),
    Input('hicool_path', 'value'))
def update_cells(hicool_path):
    global hc
    hc = HiCool(hicool_path)
    cells = sorted(list(hc.show("scool/cells")))
    field = sorted(
        list(hc.show("scool/cells/"+cells[0]+"/pixels") - set(["bin1_id", "bin2_id"])))
    features = sorted(
        list(hc.show("scool/cells/"+cells[0]+"/bins") - set(["chrom", "start", "end"])))
    chroms = hc.chroms.name
    # print(field)
    return cells, cells[0], cells, None, field, features, chroms, chroms[0]

@callback(
    Output('hicool_start1', 'value'),
    Output('hicool_end1', 'max'),
    Output('hicool_end1', 'value'),
    Output('hicool_end1', 'placeholder'),
    Input('hicool_chrom1', 'value'))
def update_cells(hicool_chrom):
    length = hc.chroms.set_index('name').loc[hicool_chrom][0]
    return 0, length, length, 'end:'+str(length)


@callback(
    Output('compare_heatmap', 'figure'),
    Output('hicool_zmax1', 'max'),
    Input('hicool_cell1', 'value'),
    Input('hicool_cell2', 'value'),
    State('hicool_chrom1', 'value'),
    State('hicool_start1', 'value'),
    State('hicool_end1', 'value'),
    Input('hicool_field1', 'value'),
    Input('hicool_features', 'options'),
    Input('hicool_zmax1', 'value'),
    Input('chrom_button', 'n_clicks'),
)
def update_graph(cell1,
                 cell2, 
                 chrom, 
                 start: int, 
                 end: int, 
                 field, 
                 features,
                 zmax, 
                 chrom_click):
    binsize = hc.info["bin-size"]
    verb = ""
    if end is None:
        end = hc.chroms.set_index('name').loc[chrom][0]
    bin_len = (int(end) - int(start)) // binsize
    if bin_len > 3000:
        end = int(start) + 3000 * binsize
        verb = "Only support matrix length less 3000 bins!\n Reshaped matrix to (3000, 3000)."
    region = f"{chrom}:{start}-{end}"
    if cell2 is not None:
        cell1_path = hc.root + "::scool/cells/" + cell1
        cell2_path = hc.root + "::scool/cells/" + cell2
        mat1 = load_cool_region(cell1_path, region, as_pixels=False, field=field)
        mat2 = load_cool_region(cell2_path, region, as_pixels=False, field=field)
        matrix = np.triu(mat1) + np.tril(mat2)
    else:
        cell_path = hc.root + "::scool/cells/" + cell1
        matrix = load_cool_region(cell_path, region, as_pixels=False, field=field)
    label = list(str(i*binsize//1000000) +
                 'M' for i in range(int(start)//binsize, int(end)//binsize))
    if matrix.shape[0] > len(label):
        label += [str((int(end)//1000000)+1) + 'M']
    if zmax is None:
        fig = px.imshow(matrix, x=label, y=label,
                        color_continuous_scale="Reds", aspect='equal')
    else:
        fig = px.imshow(matrix, x=label, y=label, zmax=zmax,
                        color_continuous_scale="Reds", aspect='equal')
    fig.update_layout(
        autosize=False,
        title=region, title_x=0.5,
        width=800, height=800,
        annotations=[dict(text=verb, 
                          showarrow=False, 
                          xref='paper', yref='paper',
                          x=0.5, y=-0.18, 
                          font=dict(size=14))],
        coloraxis_colorbar=dict(x=-0.15)
    )
    fig.update_xaxes(title_text=cell1)
    fig.update_yaxes(title_text=cell2, side='right')
    return fig, int(np.quantile(matrix,0.995))

