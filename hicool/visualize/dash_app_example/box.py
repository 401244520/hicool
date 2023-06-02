import dash
import dash_core_components as dcc
import dash_html_components as html

app = dash.Dash(__name__)

app.layout = html.Div(children=[
    html.H1(children='运算函数'),

    html.Div(children='''
        请选择要运行的选项：
    '''),

    dcc.Dropdown(
        options=[
            {'label': '选项1', 'value': 'option1'},
            {'label': '选项2', 'value': 'option2'},
            {'label': '选项3', 'value': 'option3'}
        ],
        value='option1'
    ),

    html.Button('运行', id='button'),

    html.Div(id='output')
])

@app.callback(
    dash.dependencies.Output('output', 'children'),
    [dash.dependencies.Input('button', 'n_clicks')],
    [dash.dependencies.State('dropdown', 'value')])
def run_function(n_clicks, value):
    if n_clicks is not None:
        if value == 'option1':
            # 运行选项1的函数
            return '选项1的结果'
        elif value == 'option2':
            # 运行选项2的函数
            return '选项2的结果'
        elif value == 'option3':
            # 运行选项3的函数
            return '选项3的结果'

if __name__ == '__main__':
    app.run_server(debug=True)
