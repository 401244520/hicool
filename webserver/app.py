from dash import Dash, dcc, html, Input, Output
from layout_overview import hicool_view
from layout_process import process_ui
from layout_heatmap import hicool_matrix
from layout_bigwig import bw_track
from layout_diff import hicool_diff
# from layout_feature import hicool_features
# from layout_graph import multi_demo

app = Dash(__name__, assets_folder='./assets',
           suppress_callback_exceptions=True)
datasets = "../data"


app.layout = html.Div(className='layui-layout layui-layout-admin', children=[
    html.Div(className='layui-header', children=[
        html.H1('HiCool', className='layui-logo  layui-bg-black'),
        html.Ul(className='layui-nav layui-layout-left', children=[
            html.Li([
                dcc.Link('Home', href='/home',
                         className='layui-icon layui-icon-home'),
            ], className='layui-nav-item'),
            html.Li([
                dcc.Link('View', href='/view',
                         className='layui-icon layui-icon-list'),
            ], className='layui-nav-item'),
            html.Li([
                dcc.Link('Search', href='/search',
                         className='layui-icon layui-icon-search'),
            ], className='layui-nav-item'),
            html.Li([
                dcc.Link('Document', href='/document',
                         className='layui-icon layui-icon-file'),
            ], className='layui-nav-item'),
        ]),
    ]),
    html.Div([
        html.Div([
            html.Nav([
                dcc.Link(' Overview', href='/tab-view',
                         className='layui-icon layui-icon-home'),
                dcc.Link(' Function', href='/tab-func',
                         className='layui-icon layui-icon-code-circle'),
                dcc.Link(' Heatmap', href='/tab-heatmap',
                         className='layui-icon layui-icon-picture'),
                dcc.Link(' BigWig', href='/tab-bigwig',
                         className='layui-icon layui-icon-chart-screen'),
                dcc.Link(' Feature', href='/tab-feature',
                         className='layui-icon layui-icon-template-1'),
                dcc.Link(' Differential', href='/tab-diff',
                         className='layui-icon layui-icon-template'),
                dcc.Link(' Modeling', href='/tab-multi',
                         className='layui-icon layui-icon-share'),
            ], className='layui-nav-item layui-nav-itemed')
        ], className='layui-nav layui-nav-tree '),
    ], className='layui-side layui-bg-black'),

    html.Div([
        dcc.Location(id='url', refresh=False),
        html.Div(id='page-content', className="layui-elem-quote layui-text")
    ], className='layui-body'),

    html.Div([
        html.P('版权所有 © 2023 数据库', className='layui-footer-self')
    ], className='layui-footer'),
])


@app.callback(Output('page-content', 'children'),
              Input('url', 'pathname'))
def render_page_content(pathname):
    if pathname == '/tab-heatmap':
        return hicool_matrix(datasets)
    elif pathname == '/tab-bigwig':
        return bw_track(datasets)
    elif pathname == '/tab-func':
        return process_ui()
    elif pathname == '/tab-feature':
        return #hicool_features(datasets)
    elif pathname == '/tab-view':
        return hicool_view(datasets)
    elif pathname == '/tab-diff':
        return hicool_diff(datasets)
    elif pathname == '/tab-multi':
        return #multi_demo()
    else:
        return html.Div([
            html.H1('Welcome!')
        ])


if __name__ == '__main__':
    app.run_server(debug=True, port=56666, dev_tools_props_check=False)
