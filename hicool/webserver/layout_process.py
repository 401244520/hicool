from dash import dcc, html, Input, Output, State,callback,ALL
import inspect
import textwrap

try:
    from hicool.tools import HiCool
    from hicool.function.conformation import compartment_decomposition,tad_insulation
except:
    import sys
    sys.path.append("/home/wzl/Workspace/HiCool")
    from hicool.tools import HiCool
    from hicool.function.conformation import compartment_decomposition,tad_insulation

func_api = {
    'compartment':compartment_decomposition,
    'tad':tad_insulation,
}
def process_ui():
   return html.Div([
        dcc.Dropdown(
            id='function-dropdown',
            options=list(func_api.keys()),
            value= 'compartment',
        ),
        html.Div(id='input-container'),
        html.Button('Run', id='run-button'),
        html.Div(id='output')
    ])

def my_function(a, b=2, c='hello'):
    pass

@callback(
    Output('input-container', 'children'),
    Input('function-dropdown', 'value')
)
def update_input_container(func_name):
    func = func_api[func_name]
    sig = inspect.signature(func)
    docstring = f'''```{func.__doc__} ```'''
    doc = [dcc.Markdown(docstring)]
    inp = []
    for name, param in sig.parameters.items():
        if param.default is param.empty:
            default_value = ''
        else:
            default_value = param.default

        inp.append(html.Div([
            html.Div([
                html.Label(name),
            ],style={'display': 'inline-block', 'width': '40%'}),
            html.Div([
                dcc.Input(id=name, value=default_value)
            ],style={'display': 'inline-block', 'width': '40%','textAlign': 'center'}),
        ]))

    return html.Div([
            html.Div(inp,style={'display': 'inline-block', 'width': '20%'}),
            html.Div(doc,style={'display': 'inline-block'}),
        ])

@callback(
    Output('output', 'children'),
    Input('function-dropdown', 'value'),
    Input('run-button', 'n_clicks'),
    State({'index': ALL}, 'value')
)
def run_function(func_name,n_clicks, *args):
    print(args)
    if n_clicks is None:
        return ''
    else:
        func = func_api[func_name]
        sig = inspect.signature(func)
        param_names = list(sig.parameters.keys())
        params = {name: value for name, value in zip(param_names, args)}
        print(args,param_names)
        result = func(**params)
        return result