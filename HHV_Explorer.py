# -*- coding: utf-8 -*-
import dash
import re
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import pandas as pd
from dash.dependencies import Input, Output, State
import matplotlib, random
from flask_caching import Cache
import json
import uuid
from dash.exceptions import PreventUpdate
import os
import dash_bio as dashbio
import dash_table
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.tools as tls
import math

dirpath = './' 

external_stylesheets = [dbc.themes.BOOTSTRAP,'https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.config.suppress_callback_exceptions = True

cache = Cache(app.server, config={
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': 'cache-directory'
})


TIMEOUT = 1

colrange = pd.read_pickle(dirpath+"GLA_25506/GLA_25506|1000|0.pickle")
colsam = {}

for keys in [*colrange]:
    col = 'crimson'
    xcol = 0
    while col == 'crimson':
        col = colrange[keys].color.iloc[xcol]
        xcol += 1
        colsam[keys] = col

colrange = pd.read_pickle(dirpath+"DER512/DER512|1000|0.pickle")

for keys in [*colrange]:
    col = 'crimson'
    xcol = 0
    while col == 'crimson':
        col = colrange[keys].color.iloc[xcol]
        xcol += 1
        colsam[keys] = col

colrange = None

samcon = {"iciHG01277":"iciHG01277 6A - 19q",
"LF-3A":"LF-3A 6A - 19q",
"LF-1A":"LF-1A 6A - 19q",
"LEI_1501":"LEI_1501 6A - 19q",
"506-035":"506-035 6A - 19q",
"GLA_25506":"GLA_25506 6A - 19q",
"LF-2A":"LF-2A 6A - 19q",
"JHPT-G1":"JHPT-G1 6A - 19q",
"HP73F12":"HP73F12 6A - 18q",
"JHPT-D12":"JHPT-D12 6A - 18q",
"HP94B11":"HP94B11 6A - 18q",
"HP23A7":"HP23A7 6A - 18q",
"107-040":"107-040 6A - 18q",
"303-046":"303-046 6A - 17p",
"303-035":"303-035 6A - 17p",
"GLA_15137":"GLA_15137 6A - 17p",
"7A-17p13.3":"7A-17p13.3 6A - 17p",
"103-091":"103-091 6A - 17p",
"HGDP00628":"HGDP00628 6A - 17p",
"HP88D9":"HP88D9 6A - 17p",
"HP73C5":"HP73C5 6A - 17p",
"HP15A11":"HP15A11 6A - 17p",
"HP96H12":"HP96H12 6A - 17p",
"HP104A5":"HP104A5 6A - 17p",
"LF-5A":"LF-5A 6A - 17p" ,
"NC_001664.4_U1102":"NC_001664.4_U1102 6A - NI",
"SIE":"SIE 6A - NI",
"KP257584.1_AJ":"KP257584.1_AJ 6A - NI",
"CO4":"CO4 6A - NI",
"CO7":"CO7 6A - NI",
"DA":"DA 6A - NI",
"KJ123690.1_GS":"KJ123690.1_GS 6A - NI",
"CO2":"CO2 6A - NI",
"CO1":"CO1 6A - NI",
"CO3":"CO3 6A - NI",
"NA18999":"NA18999 6A - 22q",
"LEI_ALD":"LEI_ALD 6B - 11p",
"4B-11p15.5":"4B-11p15.5 6B - 11p",
"NA07022":"NA07022 6B - 11p ",
"CUM082":"CUM082 6B -  9q",
"HP12G6":"HP12G6 6B -  9q",
"GLA_35629":"GLA_35629 6B -  9q",
"GLA_34108":"GLA_34108 6B -  9q",
"GLA_29221":"GLA_29221 6B -  9q",
"HP34B2":"HP34B2 6B -  9q",
"2B-9q34.3":"2B-9q34.3 6B -  9q",
"GLA_3986":"GLA_3986 6B -  9q",
"d37":"d37 6B -  9q",
"BAN519":"BAN519 6B -  9q",
"NA10863":"NA10863 6B -  9q",
"HP58A9":"HP58A9 6B -  19q",
"iciHG00245":"iciHG00245 6B -  19q",
"COR264":"COR264 6B -  19q",
"HP24D3":"HP24D3 6B -  19q",
"HP1F1":"HP1F1 6B -  19q",
"HP12H12":"HP12H12 6B -  19q",
"1B-iciHHV-6B":"1B-iciHHV-6B 6B - 17p(minor)",
"704-021":"704-021 6B - 17p(minor)",
"704-016":"704-016 6B - 17p(minor)",
"HGDP01065":"HGDP01065 6B - 17p(major)",
"HP46B12":"HP46B12 6B - 17p(major)",
"DER512":"DER512 6B - 17p(major)",
"HGDP01077":"HGDP01077 6B - 17p(major)",
"801-018":"801-018 6B - 17p(major)",
"ORCA3835":"ORCA3835 6B - 17p(major)",
"ORCA1622":"ORCA1622 6B - 17p(major)",
"HP91B10":"HP91B10 6B - 17p(major)",
"HP63F10":"HP63F10 6B - 17p(major)",
"Z29":"Z29 6B - NI",
"NAK":"NAK 6B - NI",
"NY-233":"NY-233 6B - NI",
"NY-393":"NY-393 6B - NI",
"japan-b3":"japan-b3 6B - NI",
"japan-b7":"japan-b7 6B - NI",
"HST":"HST 6B - NI",
"japan-b2":"japan-b2 6B - NI",
"japan-a7":"japan-a7 6B - NI",
"NY-399":"NY-399 6B - NI",
"NY-350":"NY-350 6B - NI",
"NY-32":"NY-32 6B - NI",
"MAR":"MAR 6B - NI",
"NY-241":"NY-241 6B - NI",
"KYO":"KYO 6B - NI",
"ENO":"ENO 6B - NI",
"NY-338":"NY-338 6B - NI",
"NY-104":"NY-104 6B - NI"}


refs = [['U1102 (6A)','NC_001664.4_U1102'],
    ['NA18999 (6A)','NA18999'],
    ['GLA_25506 (6A)','GLA_25506'],
    ['GLA_15137 (6A)','GLA_15137'],
    ["3A-10q26 (6A)","3A-10q26"],
    ["HGDP00092 (6B)","HGDP00092"],
    ["4B-11P15-5 (6B)","4B-11p15.5"],
    ["BAN519 (6B)","BAN519"],
    ["COR264 (6B)","COR264"],
    ["DER512 (6B)","DER512"],
    ["Z29 (6B)","Z29"],
    ["1B-HHV-6B (6B)","1B-iciHHV-6B"],
    ["HGDP00813 (6B)","HGDP00813"],
    ["HST (6B)","HST"]]

app.layout = html.Div(children=[
        html.H1(children='HHV-6 Explorer',style={'text-align':'center'}),
        html.Div(html.P(children='HHV-6 explorer is a tool to compare variation across open reading frames and genomes of inherited and non-inherited HHV-6.',style={'text-align':'center'}),style={'margin':'auto','width':'1000px'}),

    html.Div([
    html.Div([
    html.Div([
        html.Div([
        html.P(children='Select a reference HHV6 genome and sliding count window size and overlap.',style={'text-align':'center'}),
        dcc.Dropdown(
            id='reference-selector',
            options=[{'label':i,'value':j} for i,j in refs],
            value=refs[1][1],
            clearable=False,
            optionHeight=25,
            style={'text-align':'center'}
        ),
        dcc.Dropdown(
            id='tile-selector',
            clearable=False,
            optionHeight=25,
            style={'text-align':'center'}

        )], style={'margin':'auto','padding':'10px','border-radius':'5px', 'border':'1px solid lightgrey'}, className="six columns"),
        html.Div(id="cladepic", className = "six columns")
    ],className="row", style={'margin-left':'5%', 'margin-right':'5%'}),
        html.Br(),




        ]),
        html.Div(html.P(children="From the dropdown menus below, choose the type or types of variation to explore. When multiple types of variation are selected, these are cumulative on the bar charts shown. Samples can be selected by integration site or added from the dropdown menu (which can be filtered by typing a known sample).  All variation counts are relative to the reference genome chosen above. The upper and lower graphs allow for easier comparison of variation between samples or clades.",style={'text-align':'center'}),style={'margin':'auto','width':'1200px'}),
        html.Div([
        html.Div([
            html.H4('Upper Graph', style={'text-align':'center'}),
                html.H5("Filter by variation type:", className="control_label"),

                dcc.Dropdown(
            id='seqvar-typeA',
            multi=True,
            optionHeight=25,
        value=['var_substitutions']
        ),
        html.H5("Filter by Sample:", className="control_label"),
        dcc.RadioItems(
            id="sample_selector",
            options=[
                {"label": "Custom ", "value": "custom"},
                {"label": "None ", "value": "none"},
                {"label": "All", "value": "all"}
            ],
            value="custom",
            labelStyle={"display": "inline-block"},
            inputStyle={"margin": "5px"},
        ),
        dcc.Dropdown(
            id='sample-list',
            multi=True,
            optionHeight=25
        ),
        html.Div(id='buttonsA')            
        ], className="six columns"),

        html.Div([
            html.H4('Lower Graph', style={'text-align':'center'}),
                    html.H5("Filter by variation type:", className="control_label"),

                    dcc.Dropdown(
            id='seqvar-typeB',
            
        multi=True,
        optionHeight=25,
        value=['var_substitutions']
        ),
        html.H5("Filter by Sample:", className="control_label"),
        dcc.RadioItems(
            id="sample_selectorB",
            options=[
                {"label": "Custom ", "value": "custom"},
                {"label": "None ", "value": "none"},
                {"label": "All", "value": "all"}
            ],
            value="custom",
            inputStyle={"margin": "5px"},
            labelStyle={"display": "inline-block"}
        ),
        dcc.Dropdown(
            id='sample-listB',
            multi=True,
            optionHeight=25
        ),
        html.Div(id='buttonsB')
            
        ], className="six columns"),
    ], className="row", style={'margin-left':'5%', 'margin-right':'5%'}),
        html.Div(id='test',children = [dcc.Loading(id='loading-seq',children=(
        html.Div(html.H4("Genome View"), style={'text-align':'center'}),
        dcc.Graph(id='testgo',style={'height':800})           
    ))]),

    html.Div([
        ])], style={'columnCount':1}),
    html.Div(id='genes',children=[
        html.Div(html.H4("Gene View"), style={'text-align':'center'}),
        dcc.Dropdown(
            id='genelist',
            value='DR1',
            searchable=False, style={'width':'50%', 'align':'center', 'margin':'auto', 'text-align':'center'}
        ),
        html.Div(id='test2',children = [dcc.Loading(id='loading-gene',children=(
        dcc.Graph(id='testgo2',style={'height':600})           
    ))]),
        html.Div(id='alignment',children=[
		html.Div(html.H4("Gene Alignment"), style={'text-align':'center'}),
		html.Div(html.P(children='The alignment can be scrolled by clicking and dragging within the highlighted area in the full view at the bottom, it can also be resized by clicking and dragging on the white handles of the highlighted area.', style={'text-align':'center'}),style={'margin':'auto','width':'1000px'}),
        dashbio.AlignmentChart(
            id='alignment-chart',
            data = '>test',
            showconservation=False,
            showgap=False,
            showconsensus=False
        ) ], style={'margin-left':'5%', 'margin-right':'10%'})   
    ]
    ),
])


@app.callback(
    Output("modal", "is_open"),
    [Input("open", "n_clicks"), Input("close", "n_clicks")],
    [State("modal", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open

@cache.memoize(timeout=TIMEOUT)  # in seconds
def create_dual_fig(samA, samB,typeA,typeB,ref, tile, gene=None):
    print("REF IS:",ref)
    if gene != None:
        fig = make_subplots(rows = 2, cols = 1, specs=[[{}],[{}]],shared_xaxes=True,vertical_spacing=0.050, row_heights=[0.5,0.5])
        gpick = gene.split('.')[0]+'|'+tile.split('|')[1]+'|'+tile.split('|')[2]
        godf = pd.read_pickle(dirpath+ref+'/genes/'+gpick+'.pickle')
        total = pd.DataFrame(columns=['gaatotal'])
        colors = pd.DataFrame(columns=['c'])
    else:
        fig = make_subplots(rows = 3, cols = 1, specs=[[{}],[{}],[{}]],shared_xaxes=True,vertical_spacing=0.050, row_heights=[0.4,0.4,0.2])
        godf = pd.read_pickle(dirpath+ref+'/'+tile+'.pickle')
        total = pd.DataFrame(columns=['total'])
        colors = pd.DataFrame(columns=['c'])

    if samA == None or len(samA) == 0:
        samA = [next(iter(godf))]
        while samA[0] == ref:
            samA = [next(iter(godf))]
    if samB == None or len(samB) ==  0:
        samB = [next(iter(godf))]

    for x in samA:
        if x not in list(godf.keys()):
            samA.remove(x)
    for x in samB:
        if x not in list(godf.keys()):
            samB.remove(x)

    if samA == None or len(samA) == 0:
        samA = [next(iter(godf))]
    if samB == None or len(samB) ==  0:
        samB = [next(iter(godf))]
    max = 0 # to allow both graphs to have same range for y
    for sample in samA:
        firstatt = True
        for att in typeA:
            if firstatt:
                total[sample] = godf[sample][att]
                colors[sample] = godf[sample]['color']
            else: 
                total[sample] = total[sample] + godf[sample][att]
            firstatt = False
        colors[sample] = colors[sample].where(colors[sample]=='crimson',colsam[sample])
        total[sample] = total[sample].where((total[sample]!=0) | (godf[sample]['color']!='crimson'),-1)
        max = total[sample].max() if total[sample].max() > max else max
        fig.append_trace({'x':godf[sample].index.values, 'y':total[sample].values, 'marker_color':colors[sample].values,'marker_line_width':0,'type':'bar', 'name':sample},1,1)
    for sample in samB:
        firstatt = True
        for att in typeB:
            if firstatt:
                total[sample] = godf[sample][att]
                colors[sample] = godf[sample]['color']
            else: 
                total[sample] = total[sample] + godf[sample][att]
            firstatt = False
        colors[sample] = colors[sample].where(colors[sample]=='crimson',colsam[sample])
        total[sample] = total[sample].where((total[sample]!=0) | (godf[sample]['color']!='crimson'),-1)
        max = total[sample].max() if total[sample].max() > max else max
        fig.append_trace({'x':godf[sample].index.values, 'y':total[sample].values, 'marker_color':colors[sample].values,'marker_line_width':0,'type':'bar', 'name':sample},2,1)
    
    if gene == None:
        gd = pd.read_pickle(dirpath+ref+'/genes/genes.direction')
        gl = list(gd.keys())
        tpos = 'bottom right'
        spacer = 0
        for g in gl:
            if tpos == 'bottom right':
                tpos = 'top right'
                spacer =  0
            else:
                spacer =  5
                tpos = 'bottom right'
            godf2 = pd.read_pickle(dirpath+ref+'/genes/'+g+'|500|0.pickle')
            gaatotal = pd.DataFrame(columns=['gaatotal'])
            samp = [ref]
            for sample in samp:
                firstatt = True
                for att in typeA:
                    if firstatt:
                        gaatotal[sample] = ((godf2[sample][att]+1)/(godf2[sample][att]+1))+spacer
                    else:
                        gaatotal[sample] = ((godf2[sample][att]+1)/(godf2[sample][att]+1))+spacer
                    firstatt = False
                gaatotal[sample] = gaatotal[sample].where((gaatotal[sample]!=0) | (godf2[sample]['color']!='crimson'),-1)
                fig.append_trace(go.Scatter(x=godf2[sample].index.values, y=gaatotal[sample].values, showlegend=False, text=[g],mode='lines+markers+text',hoverinfo='none', textposition=tpos,textfont_size=10, line_width = 5,marker_color='black'),3,1)
    if gene == None:
        fig['layout'].update( title='Cumulative variation for selected samples across the HHV-6 genome', showlegend=False, barmode='group', bargap=0, plot_bgcolor='#fff', yaxis_fixedrange=True,yaxis2_fixedrange=True, yaxis3_fixedrange=True)
        fig.update_yaxes(title_text="Genes",visible=False, row=3, col=1)
        fig.update_xaxes(title_text="BP", row=3, col=1)
    else:
        fig['layout'].update( title='Cumulative variation for selected samples across the selected gene', showlegend=False, barmode='group', bargap=0, plot_bgcolor='#fff', yaxis_fixedrange=True,yaxis2_fixedrange=True)
        fig.update_xaxes(title_text="BP", row=2, col=1)

    max = 5 * round((max+2.4)/5)
    fig.update_yaxes(title_text="Count", range=[-1,max], row=1, col=1)
    fig.update_yaxes(title_text="Count", range=[-1,max], row=2, col=1)
    fig.update_yaxes(gridcolor='lightgrey', zerolinecolor='black')
    return fig

@app.callback(
    [Output('buttonsA','children'),
    Output('buttonsB','children')],
    [Input('sample-list', 'value'),
    Input('sample-listB', 'value'),
    Input('reference-selector', 'value'),
    Input('tile-selector','value')]
)
def butcols(samA,samB, ref, tile):
    if samA == None or samB == None or  ref == None or tile == None or not os.path.exists(dirpath+ref+'/'+tile+'.pickle'):                         
        raise PreventUpdate 

    chilA = []
    chilB = []
   
    for keys in samA:
        chilA.append(html.Button(samcon[keys], style={'text-transform': 'none','background-color':colsam[keys],'font-size':'15px','color':'#ffffff'}))
    for keys in samB:
        chilB.append(html.Button(samcon[keys], style={'text-transform': 'none','background-color':colsam[keys],'font-size':'15px','color':'#ffffff'}))  

    return chilA,chilB

@app.callback(
    [Output('testgo','figure')],
    [
    Input('sample-list', 'value'),
    Input('sample-listB', 'value'),
    Input('seqvar-typeA','value'),
    Input('seqvar-typeB','value'),
    Input('reference-selector', 'value'),
    Input('tile-selector','value')
    ]
)
def updatedualgraph(samA, samB,typeA,typeB,ref, tile):
    if tile == None or tile == None or not os.path.exists(dirpath+ref+'/'+tile+'.pickle'):                         
        raise PreventUpdate 

    fig = create_dual_fig(samA, samB,typeA,typeB,ref, tile)

    return [fig]

@app.callback(
    [Output('testgo2','figure')],
    [
    Input('sample-list', 'value'),
    Input('sample-listB', 'value'),
    Input('seqvar-typeA','value'),
    Input('seqvar-typeB','value'),
    Input('reference-selector', 'value'),
    Input('tile-selector','value'),
    Input('genelist','value')
    ]
)
def updatedualgraph2(samA, samB,typeA,typeB,ref, tile, gene):
    if ref == None or tile == None or not os.path.exists(dirpath+ref+'/'+tile+'.pickle'):                         
        raise PreventUpdate 
    fig = create_dual_fig(samA, samB,typeA,typeB,ref, tile,gene)
    return [fig]

@cache.memoize(timeout=TIMEOUT)  # in seconds
def process_dropdown_lists(ref,tile):
    gd = pd.read_pickle(dirpath+ref+'/genes/genes.direction')
    gl = list(gd.keys())
    def atoi(text):
        return int(text) if text.isdigit() else text
    def natural_keys(text):
        return [ atoi(c) for c in re.split('(\d+)', text)]
    gl.sort(key=natural_keys)
    gsopt = [{'label':i,'value':i} for i in gl]
    return [gsopt,gsopt[0]['value']]

#load files based on selected reference
@app.callback(
    [
    Output('seqvar-typeA','options'),
    Output('seqvar-typeA','value'),
    Output('seqvar-typeB','options'),
    Output('seqvar-typeB','value'),
    Output('sample-list','options'),
    Output('sample-listB','options'),
    Output('genelist','options'),
    Output('genelist','value'),
    Output('sample_selector','options'),
    Output('sample_selector','value'),
    Output('sample_selectorB','options'),
    Output('sample_selectorB','value'),
    Output('cladepic','children')],
    [Input('reference-selector', 'value'),
    Input('tile-selector','value'),
    State('sample_selector','value'),
    State('sample_selectorB','value')]
)
def read_data(ref,tile,sa,sb):
    print('changing Reference or Window',ref,tile)
    if ref == None or tile == None or not os.path.exists(dirpath+ref+'/'+tile+'.pickle'):
        raise PreventUpdate
    else:
        cpick = 'HHV-6A_network_forexplorer.png'
        ctitle = "HHV-6A clades"
        dlists = process_dropdown_lists(ref,tile)
        if ref in ['NA18999','NC_001664.4_U1102','107-040','GLA_25506','GLA_15137','3A-10q26']:
            ctitle = 'HHV-6A clades'
            cpick = 'HHV-6A_network_forexplorer.png'
            if sa not in ['custom','19qA','18q','17p','nonA']:
                sa = 'custom'
            if sb not in ['custom','19qA','18q','17p','nonA']:
                sb = 'custom'
            ssopt = [
                    {"label": "Custom ", "value": "custom"},
                    {"label": "19q ", "value": "19qA"},
                    {"label": "18q ", "value": "18q"},
                    {"label": "17p ", "value": "17p"},
                    {"label": "Non-inherited", "value": "nonA"}
                ]
            samples = [
                    {"label": "iciHG01277 6A - 19q ", "value": "iciHG01277"},
                    {"label": "LF-3A 6A - 19q ", "value": "LF-3A"},
                    {"label": "LF-1A 6A - 19q ", "value": "LF-1A"},
                    {"label": "LEI_1501 6A - 19q ", "value": "LEI_1501"},
                    {"label": "506-035 6A - 19q ", "value": "506-035"},
                    {"label": "GLA_25506 6A - 19q ", "value": "GLA_25506"},
                    {"label": "LF-2A 6A - 19q ", "value": "LF-2A"},
                    {"label": "JHPT-G1 6A - 19q ", "value": "JHPT-G1"},
                    {"label": "HP73F12 6A - 18q ", "value": "HP73F12"},
                    {"label": "JHPT-D12 6A - 18q ", "value": "JHPT-D12"},
                    {"label": "HP94B11 6A - 18q ", "value": "HP94B11"},
                    {"label": "HP23A7 6A - 18q ", "value": "HP23A7"},
                    {"label": "107-040 6A - 18q ", "value": "107-040"},
                    {"label": "303-046 6A - 17p ","value": "303-046"},
                    {"label": "303-035 6A - 17p ","value": "303-035"},
                    {"label": "GLA_15137 6A - 17p ","value": "GLA_15137"},
                    {"label": "7A-17p13.3 6A - 17p ","value": "7A-17p13.3"},
                    {"label": "103-091 6A - 17p ","value": "103-091"},
                    {"label": "HGDP00628 6A - 17p ","value": "HGDP00628"},
                    {"label": "HP88D9 6A - 17p ","value": "HP88D9"},
                    {"label": "HP73C5 6A - 17p ","value": "HP73C5"},
                    {"label": "HP15A11 6A - 17p ","value": "HP15A11"},
                    {"label": "HP96H12 6A - 17p ","value": "HP96H12"},
                    {"label": "HP104A5 6A - 17p ","value": "HP104A5"},
                    {"label": "LF-5A 6A - 17p ","value": "LF-5A"},
                    {"label": "NC_001664.4_U1102 6A - NI ", "value": "NC_001664.4_U1102"},
                    {"label": "SIE 6A - NI ", "value": "SIE"},
                    {"label": "KP257584.1_AJ 6A - NI ", "value": "KP257584.1_AJ"},
                    {"label": "CO4 6A - NI ", "value": "CO4"},
                    {"label": "CO7 6A - NI ", "value": "CO7"},
                    {"label": "DA 6A - NI ", "value": "DA"},
                    {"label": "KJ123690.1_GS 6A - NI ", "value": "KJ123690.1_GS"},
                    {"label": "CO2 6A - NI ", "value": "CO2"},
                    {"label": "CO1 6A - NI ", "value": "CO1"},
                    {"label": "CO3 6A - NI ", "value": "CO3"},
                    {"label": "NA18999 6A - 22q ", "value": "NA18999"}                
            ]
        else:
            cpick = 'HHV-6B_network_forexplorer.png'
            ctitle = 'HHV-6B clades'
            if sa not in ['custom','11p','9q','19qB','17pm','17pmj','nonB']:
                sa = 'custom'
            if sb not in ['custom','11p','9q','19qB','17pm','17pmj','nonB']:
                sb = 'custom'
            ssopt = [
                    {"label": "Custom ", "value": "custom"},
                    {"label": "11p ", "value": "11p"},
                    {"label": "9q ", "value": "9q"},
                    {"label": "19q ", "value": "19qB"},
                    {"label": "17p(minor) ", "value": "17pm"},
                    {"label": "17p(major) ", "value": "17pmj"},
                    {"label": "Non-inherited", "value": "nonB"}
                ]
            samples = [
                    {"label": "LEI_ALD 6B - 11p ","value": "LEI_ALD"},
                    {"label": "4B-11p15.5 6B - 11p ","value": "4B-11p15.5"},
                    {"label": "NA07022 6B - 11p  ","value": "NA07022"},
                    {"label": "CUM082 6B -  9q ","value": "CUM082"},
                    {"label": "HP12G6 6B -  9q ","value": "HP12G6"},
                    {"label": "GLA_35629 6B -  9q ","value": "GLA_35629"},
                    {"label": "GLA_34108 6B -  9q ","value": "GLA_34108"},
                    {"label": "GLA_29221 6B -  9q ","value": "GLA_29221"},
                    {"label": "HP34B2 6B -  9q ","value": "HP34B2"},
                    {"label": "2B-9q34.3 6B -  9q ","value": "2B-9q34.3"},
                    {"label": "GLA_3986 6B -  9q ","value": "GLA_3986"},
                    {"label": "d37 6B -  9q ","value": "d37"},
                    {"label": "BAN519 6B -  9q ","value": "BAN519"},
                    {"label": "NA10863 6B -  9q ","value": "NA10863"},
                    {"label": "HP58A9 6B -  19q ","value": "HP58A9"},
                    {"label": "iciHG00245 6B -  19q ","value": "iciHG00245"},
                    {"label": "COR264 6B -  19q ","value": "COR264"},
                    {"label": "HP24D3 6B -  19q ","value": "HP24D3"},
                    {"label": "HP1F1 6B -  19q ","value": "HP1F1"},
                    {"label": "HP12H12 6B -  19q ","value": "HP12H12"},
                    {"label": "1B-iciHHV-6B 6B - 17p(minor) ","value": "1B-iciHHV-6B"},
                    {"label": "704-021 6B - 17p(minor) ","value": "704-021"},
                    {"label": "704-016 6B - 17p(minor) ","value": "704-016"},
                    {"label": "HGDP01065 6B - 17p(major) ","value": "HGDP01065"},
                    {"label": "HP46B12 6B - 17p(major) ","value": "HP46B12"},
                    {"label": "DER512 6B - 17p(major) ","value": "DER512"},
                    {"label": "HGDP01077 6B - 17p(major) ","value": "HGDP01077"},
                    {"label": "801-018 6B - 17p(major) ","value": "801-018"},
                    {"label": "ORCA3835 6B - 17p(major) ","value": "ORCA3835"},
                    {"label": "ORCA1622 6B - 17p(major) ","value": "ORCA1622"},
                    {"label": "HP91B10 6B - 17p(major) ","value": "HP91B10"},
                    {"label": "HP63F10 6B - 17p(major) ","value": "HP63F10"},
                    {"label": "Z29 6B - NI ","value": "Z29"},
                    {"label": "NAK 6B - NI ","value": "NAK"},
                    {"label": "NY-233 6B - NI ","value": "NY-233"},
                    {"label": "NY-393 6B - NI ","value": "NY-393"},
                    {"label": "japan-b3 6B - NI ","value": "japan-b3"},
                    {"label": "japan-b7 6B - NI ","value": "japan-b7"},
                    {"label": "HST 6B - NI ","value": "HST"},
                    {"label": "japan-b2 6B - NI ","value": "japan-b2"},
                    {"label": "japan-a7 6B - NI ","value": "japan-a7"},
                    {"label": "NY-399 6B - NI ","value": "NY-399"},
                    {"label": "NY-350 6B - NI ","value": "NY-350"},
                    {"label": "NY-32 6B - NI ","value": "NY-32"},
                    {"label": "MAR 6B - NI ","value": "MAR"},
                    {"label": "NY-241 6B - NI ","value": "NY-241"},
                    {"label": "KYO 6B - NI ","value": "KYO"},
                    {"label": "ENO 6B - NI ","value": "ENO"},
                    {"label": "NY-338 6B - NI ","value": "NY-338"},
                    {"label": "NY-104 6B - NI ","value": "NY-104"}
            ]
        
        variation = [
            {'label': 'substitutions', 'value': 'var_substitutionsAll'},  
            {'label': 'insertions', 'value': 'var_insertions'}, 
            {'label': 'deletions', 'value': 'var_deletions'}, 
            {'label': 'coding substitutions', 'value': 'var_coding_substitutionsAll'}, 
            {'label': 'coding insertions', 'value': 'var_coding_insertions'}, 
            {'label': 'coding deletions', 'value': 'var_coding_deletions'}, 
            {'label': 'missense', 'value': 'aa_missenseAll'},  
            {'label': 'nonsense', 'value': 'aa_nonsense'}, 
            {'label': 'truncation', 'value': 'aa_truncation'}, 
            {'label': 'AA insertion', 'value': 'aa_gene-insertion'}, 
            {'label': 'AA deletions', 'value': 'aa_deletions'}, 
            {'label': 'frameshift', 'value': 'aa_frameshift'}, 
            {'label': 'frameshift-missense', 'value': 'aa_frameshift-missense'}, 
            {'label': 'splice-site', 'value': 'aa_splice-site'}, 
            {'label': 'no cds', 'value': 'aa_nocds'}, 
            {'label': 'run-on', 'value': 'aa_run-on'}
        ]
        modal = html.Div(
            [
                dbc.Button("Enlarge Image", id="open"),
                dbc.Modal(
                    [
                        dbc.ModalHeader(ctitle),
                        dbc.ModalBody(html.Img(src=app.get_asset_url(cpick), style={'height':900})),
                        dbc.ModalFooter(
                            dbc.Button("Close", id="close", className="ml-auto")
                        ),
                    ],
                    id="modal",
                    size="xl"
                ),
            ]
        )
        cp = [html.Img(src=app.get_asset_url(cpick), style={'height':400}),modal]
        return variation,['var_substitutionsAll'],variation,['var_substitutionsAll'],samples,samples,dlists[0],dlists[1],ssopt,sa,ssopt,sb,cp

@app.callback(
    [Output('tile-selector', 'options'),
    Output('tile-selector','value')],
    [Input('reference-selector', 'value')]
)
def update_reference(dropdown):
    if dropdown == None:
        raise PreventUpdate
    else:
        countsfolder = os.listdir(dirpath+dropdown)
        counts = []
        for file in countsfolder:
            if 'pickle' in file:
                process = file.replace('.pickle','')
                process = process.split('|')
                if len(process) > 1:
                    label = 'window:'+process[1]+' overlap:'+process[2]
                    counts.append([int(process[1]),int(process[2]),label,file.replace('.pickle','')])
        counts.sort(key=lambda x: (x[0],x[1]))
        opt = [{'label':i,'value':j} for a,b,i,j in counts]
        val = opt[0]['value']
        return opt,val

@app.callback(
   [Output("sample-list", "value"), 
   Output("sample-listB", "value"), ],
    [Input("sample_selector", "value"),
    Input("sample_selectorB", "value"),
    Input("sample-list","options"),
    Input("sample-listB","options"),
    Input('reference-selector', 'value')]
)
def sample_status(selectorA,selectorB,optA,optB,ref):
    if optA == None or optB == None:
        raise PreventUpdate
    else:
        retA = []
        retB = []

        if selectorA == "custom":
            if ref in ['NA18999','NC_001664.4_U1102','107-040','GLA_25506','GLA_15137','3A-10q26']:
                retA = ["LF-3A"]
            else:
                retA = ["LEI_ALD"]
        elif selectorA == "19qA":
            retA = ['iciHG01277','LF-3A','LF-1A','LEI_1501','506-035','GLA_25506','LF-2A','JHPT-G1']
        elif selectorA == "18q":
            retA = ['HP73F12','JHPT-D12','HP94B11','HP23A7','107-040']
        elif selectorA == "17p":
            retA = ['303-046','303-035','GLA_15137','7A-17p13.3','103-091','HGDP00628','HP88D9','HP73C5','HP15A11','HP96H12','HP104A5','LF-5A']
        elif selectorA == "nonA":
            retA = ['NC_001664.4_U1102','SIE','KP257584.1_AJ','CO4','CO7','DA','KJ123690.1_GS','CO2','CO1','CO3']
        elif selectorA == "22q":
            retA = ['NA18999']
        elif selectorA == "11p":
            retA =  ['LEI_ALD','4B-11p15.5','NA07022']
        elif selectorA == "9q":
            retA =  ['CUM082','HP12G6','GLA_35629','GLA_34108','GLA_29221','HP34B2','2B-9q34.3','GLA_3986','d37','BAN519','NA10863']
        elif selectorA == "19qB":
            retA =  ['HP58A9','iciHG00245','COR264','HP24D3','HP1F1','HP12H12']
        elif selectorA == "17pm":
            retA =  ['1B-iciHHV-6B','704-021','704-016']
        elif selectorA == "17pmj":
            retA =  ['HGDP01065','HP46B12','DER512','HGDP01077','801-018','ORCA3835','ORCA1622','HP91B10','HP63F10']
        elif selectorA == "nonB":
            retA =  ['Z29','NAK','NY-233','NY-393','japan-b3','japan-b7','HST','japan-b2','japan-a7','NY-399','NY-350','NY-32','MAR','NY-241','KYO','ENO','NY-338','NY-104']

        if selectorB == "custom":
            if ref in ['NA18999','NC_001664.4_U1102','107-040','GLA_25506','GLA_15137','3A-10q26']:
                retB = ["LF-3A"]
            else:
                retB = ["LEI_ALD"]
        elif selectorB == "19qA":
            retB = ['iciHG01277','LF-3A','LF-1A','LEI_1501','506-035','GLA_25506','LF-2A','JHPT-G1']
        elif selectorB == "18q":
            retB = ['HP73F12','JHPT-D12','HP94B11','HP23A7','107-040']
        elif selectorB == "17p":
            retB = ['303-046','303-035','GLA_15137','7A-17p13.3','103-091','HGDP00628','HP88D9','HP73C5','HP15A11','HP96H12','HP104A5','LF-5A']
        elif selectorB == "nonA":
            retB = ['NC_001664.4_U1102','SIE','KP257584.1_AJ','CO4','CO7','DA','KJ123690.1_GS','CO2','CO1','CO3']
        elif selectorB == "22q":
            retB = ['NA18999']
        elif selectorB == "11p":
            retB =  ['LEI_ALD','4B-11p15.5','NA07022']
        elif selectorB == "9q":
            retB =  ['CUM082','HP12G6','GLA_35629','GLA_34108','GLA_29221','HP34B2','2B-9q34.3','GLA_3986','d37','BAN519','NA10863']
        elif selectorB == "19q":
            retB =  ['HP58A9','iciHG00245','COR264','HP24D3','HP1F1','HP12H12']
        elif selectorB == "17pm":
            retB =  ['1B-iciHHV-6B','704-021','704-016']
        elif selectorB == "17pmj":
            retB =  ['HGDP01065','HP46B12','DER512','HGDP01077','801-018','ORCA3835','ORCA1622','HP91B10','HP63F10']
        elif selectorB == "nonB":
            retB =  ['Z29','NAK','NY-233','NY-393','japan-b3','japan-b7','HST','japan-b2','japan-a7','NY-399','NY-350','NY-32','MAR','NY-241','KYO','ENO','NY-338','NY-104']

        return retA,retB


@app.callback(
    Output('alignment-chart', 'data'),
    [Input('genelist', 'value'),
    Input('reference-selector', 'value')]
)
def update_aligment(dropdown,ref):
    if dropdown == None or ref == None:
        raise PreventUpdate
    else:
        with open(dirpath+ref+'/genes/'+dropdown+'.fasta') as fastafile:
            newdata = fastafile.read()
        return newdata


if __name__ == '__main__':
    app.run_server(debug=True)

print('finished')