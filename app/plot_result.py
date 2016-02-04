__author__ = 'liwang'

import pandas as pd
import numpy as np
from plotly.offline import plot
import plotly.graph_objs as go
from plotly import tools
from plotly.tools import FigureFactory as FF
from plotly.tools import _Table as Table


def get_text(row):
    return 'Drug A: '+row['Drug A target']+'<br>'+'Drug B: '+row['Drug B target']

def plot_ready(df_prob, df_drug):

    # map drug IDs to (1) a unique integer value for plotting, (2) its target genes
    drugID2target=dict(zip(df_drug['DrugID'],df_drug['Target_genes']))

    combID_list=[ i.split('.') for i in df_prob['Combination_ID']]
    df_ready=pd.DataFrame([[i,j,drugID2target[i],drugID2target[j]] for i,j in combID_list],columns=['Drug_A','Drug_B','Drug A target','Drug B target'])

    labels=pd.DataFrame(df_ready.apply(get_text,axis=1),columns=['Label'])
    df_prob['Prob']=[round(i,4) for i in df_prob['Prob']]
    prob=pd.DataFrame(df_prob['Prob']).reset_index()

    df_ready=pd.concat([df_ready,prob,labels],axis=1)

    return df_ready

def plot_result(df_prob, df_drug,cell_name):


    df_hbar=plot_ready(df_prob.sort_values(by=['Prob']).iloc[-20:,:],df_drug)

    hbar_label=[i[0]+' & '+i[1] for i in zip(df_hbar['Drug_A'],df_hbar['Drug_B'])]
    hbar_annotation=[round(i,3) for i in df_hbar['Prob']]

    trace1 = go.Bar(x=df_hbar['Prob'],
                    y=hbar_label,
                    text=df_hbar['Label'],
                    marker=dict(color='rgba(50, 171, 96, 0.6)',
                                line=dict(color='rgba(50, 171, 96, 1.0)',width=1,),
                                ),
                    orientation='h'
                    )
    data=[trace1]


    # Update the margins to add a title and see graph x-labels.
    annotation_spacing=list(df_hbar['Prob'])[-1]*0.05
    annotation=[dict(xref='x2',
                 yref='y2',
                 x=xd+annotation_spacing,
                 y=yd,
                 text=textd,
                 font=dict(family='Arial', size=12,color='rgb(50, 171, 96)'),
                 showarrow=False
                 ) for xd,yd,textd in zip(df_hbar['Prob'],hbar_label,hbar_annotation)
            ]

    layout=dict(
            title='Predicted synergy probability for drug combinations in cell line '+cell_name +'<br>Validated synergy measures have Probability = 1.0',
            titlefont=dict(size=20),
            showlegend=False,
            height=600,
            width=1200,
            margin=go.Margin(l=280,r=280,b=80,t=100,pad=5),paper_bgcolor='rgb(243, 243, 243)',
            plot_bgcolor='rgb(243, 243, 243)',
            hovermode='y',
            xaxis=dict(title='Synergy Probability'),
            annotations = annotation
                )
    fig = go.Figure(data=data, layout=layout)

    return plot(fig, output_type = 'div', include_plotlyjs = False,show_link=False)
