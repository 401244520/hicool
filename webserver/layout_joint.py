import glob
import plotly.express as px
import pandas as pd
from tqdm import tqdm
import numpy as np
from dash import html, dcc, Input, Output, callback, State, callback_context

from hicool.tools import HiCool
from hicool.function.dataloader import load_cool_region


def hicool_matrix():
    hicool_list = glob.glob("datasets/uploads/*hicool")
    # print(str(hicool_list) + " \n")
    return html.Div([
        html.Div([
            html.Div([
                html.Div([
                    dcc.Dropdown(
                        id='hicool_path',
                        options=hicool_list,
                        value=hicool_list[0],
                    ),
                ],style={'display': 'inline-block','textAlign': 'center', 'width': '80%'}),
                html.Div([
                    html.Button(id='hicool_button',
                                children='Show All Cells', style={'height': '30px'}),
                ], style={'display': 'inline-block', 'width': '20%','position': 'relative', 'bottom': '14px'}),
            ]),
            html.Div([
                html.Div([
                    dcc.Dropdown(id='hicool_chrom')
                ], style={'display': 'inline-block', 'width': '20%', 'height': '22px'}),
                html.Div([
                    dcc.Input(id='hicool_start', placeholder='start:0',
                                type='int', min=0, value=0, style={'height': '30px'}),
                    dcc.Input(id='hicool_end', placeholder='end:',
                                type='int', style={'height': '30px'}),
                    html.Button(id='chrom_button',
                                children='Show Genome Range', style={'height': '30px'}),
                ], style={'display': 'inline-block', 'width': '80%'}),
            ]),
            html.Div([
                html.Div([
                    dcc.RadioItems(id='hicool_field', value='count'),
                ]),
                html.Div([
                    dcc.Graph(id='hicool_heatmap'),
                ], style={'display': 'inline-block', 'width': '82%'}),
                html.Div([
                    dcc.Slider(id='hicool_zmax', min=0, max=100, vertical=True,verticalHeight=650),
                ], style={ 'display': 'inline-block','height':'750px'}),
            ]),
        ],style={'display': 'inline-block', 'width': '70%'}),
        
        html.Div([
            dcc.Dropdown(id='hicool_cell', multi=True, maxHeight=500)
        ],style={'display': 'inline-block', 'width': '30%',
                 'position': 'relative', 'right': '100px','height':'905px'}),
    ])

@callback(
    Output('hicool_cell', 'options'),
    Output('hicool_cell', 'value'),
    Output('hicool_field', 'options'),
    Output('hicool_chrom', 'options'),
    Output('hicool_chrom', 'value'),
    Input('hicool_path', 'value'))
def update_cells(hicool_path):
    global hc
    hc = HiCool(hicool_path)
    cells = sorted(list(hc.show("scool/cells")))
    field = sorted(
        list(hc.show("scool/cells/"+cells[0]+"/pixels") - set(["bin1_id", "bin2_id"])))
    chroms = hc.chroms.name
    # print(field)
    return cells, [cells[0]], field, chroms, chroms[0]


@callback(
    Output('hicool_start', 'value'),
    Output('hicool_end', 'max'),
    Output('hicool_end', 'value'),
    Output('hicool_end', 'placeholder'),
    Input('hicool_chrom', 'value'))
def update_cells(hicool_chrom):
    length = hc.chroms.set_index('name').loc[hicool_chrom][0]
    return 0, length, length, 'end:'+str(length)


@callback(
    Output('hicool_heatmap', 'figure'),
    Output('hicool_zmax', 'max'),
    Input('hicool_cell', 'value'),
    State('hicool_chrom', 'value'),
    State('hicool_start', 'value'),
    State('hicool_end', 'value'),
    Input('hicool_field', 'value'),
    Input('hicool_zmax', 'value'),
    Input('hicool_button', 'n_clicks'),
    # Input('cells_button', 'n_clicks'),
    Input('chrom_button', 'n_clicks'),
)
def update_graph(cells, 
                 chrom, 
                 start: int, 
                 end: int, 
                 field, 
                 zmax, 
                 hicool_clicks,
                #  cells_click,
                 chrom_click):
    trigger = callback_context.triggered[0]['prop_id'].split('.')[0]
    if trigger == "hicool_button":
        cells_path = hc.load_cells()
    else :
        cells_path = [hc.scool_path + "/cells/" + cell for cell in cells]
    verb = ""
    binsize = hc.info["bin-size"]
    if end is None:
        end = hc.chroms.set_index('name').loc[chrom][0]
    bin_len = (int(end) - int(start)) // binsize
    if bin_len > 3000:
        end = int(start) + 3000 * binsize
        verb = "Only support matrix length less 3000 bins!\n Reshaped matrix to (3000, 3000)."
    # cell_path = hc.scool_path + "/cells/" + cell
    matrices = []
    for cell_path in tqdm(cells_path):
        region = f"{chrom}:{start}-{end}"
        matrix = load_cool_region(cell_path, region, as_pixels=False, field=field)
        matrices.append(matrix)
    label = list(str(i*binsize//1000000) +
                 'M' for i in range(int(start)//binsize, int(end)//binsize))
    if matrix.shape[0] > len(label):
        label += [str((int(end)//1000000)+1) + 'M']
    matrices = np.nan_to_num(np.array(matrices))
    print(matrices.shape)
    if zmax is None:
        fig = px.imshow(matrices, x=label, y=label, animation_frame = 0,
                        color_continuous_scale="Reds", aspect='equal')
    else:
        fig = px.imshow(matrices, x=label, y=label, zmax=zmax, animation_frame=0,
                        color_continuous_scale="Reds", aspect='equal')
    fig.update_layout(
        autosize=False,
        title=region, title_x=0.5,
        width=850, height=800,
        annotations=[dict(text=verb, 
                          showarrow=False, 
                          xref='paper', yref='paper',
                          x=0.5, y=-0.18, 
                          font=dict(size=14))]
    )

    return fig, np.quantile(matrices,0.995)
