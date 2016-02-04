__author__ = 'liwang'

import pandas as pd
import numpy as np
from plotly.offline import plot
import plotly.graph_objs as go
from plotly import tools

def get_text(row):
    return 'Drug A: '+row['Drug_A_target']+'<br>'+'Drug B: '+row['Drug_B_target']

def plot_ready(df_prob, df_drug):

    # map drug IDs to (1) a unique integer value for plotting, (2) its target genes
    drugID2num=dict(zip(df_drug['DrugID'],np.arange(1,len(df_drug['DrugID'])+1)))
    drugID2target=dict(zip(df_drug['DrugID'],df_drug['Target_genes']))

    combID_list=[ i.split('.') for i in df_prob['Combination_ID']]
    df_ready=pd.DataFrame([[i,j,drugID2target[i],drugID2target[j],drugID2num[i],drugID2num[j]] for i,j in combID_list],columns=['Drug_A','Drug_B','Drug_A_target','Drug_B_target','xval','yval'])

    labels=pd.DataFrame(df_ready.apply(get_text,axis=1),columns=['Label'])
    prob=pd.DataFrame(df_prob['Prob']).reset_index()

    df_ready=pd.concat([df_ready,prob,labels],axis=1)

    return df_ready

def plot_bubble(df_prob, df_drug,cell_name):

    df_Bubble=plot_ready(df_prob.loc[df_prob['Prob']>0.2,:],df_drug)
    df_hbar=plot_ready(df_prob.sort_values(by=['Prob']).iloc[-20:,:],df_drug)





    hbar_label=[i[0]+' & '+i[1] for i in zip(df_hbar['Drug_A'],df_hbar['Drug_B'])]
    hbar_annotation=[round(i,3) for i in df_hbar['Prob']]
    trace1 = go.Bar(x=df_hbar['Prob'],
                    y=hbar_label,
                    text=df_hbar['Label'],
                    marker=dict(color='rgba(50, 171, 96, 0.6)',
                                line=dict(color='rgba(50, 171, 96, 1.0)',width=1,),
                                ),
                    orientation='h',
                    )
    trace2 = go.Scatter(x=df_Bubble['xval'],y=df_Bubble['yval'],text=df_Bubble['Label'],mode='markers',
                    marker=dict(size=df_Bubble['Prob']*100,
                                color=df_Bubble['Prob']*100,
                                showscale=True)
                    )

    fig = tools.make_subplots(rows=1, cols=2, specs=[[{}, {}]],subplot_titles=('', ''),vertical_spacing=1)
    fig.append_trace(trace1, 1, 1)
    fig.append_trace(trace2, 1, 2)
    annotation_spacing=list(df_hbar['Prob'])[-1]*0.1
    annotation=[dict(xref='x1',
                 yref='y1',
                 x=xd+annotation_spacing,
                 y=yd,
                 text=textd,
                 font=dict(family='Arial', size=12,color='rgb(50, 171, 96)'),
                 showarrow=False
                 ) for xd,yd,textd in zip(df_hbar['Prob'],hbar_label,hbar_annotation)
            ]
    fig['layout']['xaxis1'].update(title='Synergy Probability')
    fig['layout']['xaxis2'].update(title='Drug A',
                                   titlefont=dict(size=20),
                                   autorange=True,
                                   showgrid=True,
                                   zeroline=True,
                                   showline=False,
                                   autotick=True,
                                   ticks='',
                                   gridcolor='rgb(255, 255, 255)',
                                   showticklabels=False)

    fig['layout']['yaxis1'].update(title='')
    fig['layout']['yaxis2'].update(title='Drug B',
                                   titlefont=dict(size=20),
                                   autorange=True,
                                   showgrid=True,
                                   zeroline=True,
                                   showline=False,
                                   autotick=True,
                                   ticks='',
                                   gridcolor='rgb(255, 255, 255)',
                                   showticklabels=False)
    fig['layout'].update(title='Predicted synergy probability for drug combinations in cell line '+cell_name +'<br>Validated synergy measures have Probability = 1.0',
                         titlefont=dict(size=20),
                         showlegend=False,
                         height=600,
                         width=1200,
                         margin=go.Margin(l=250,r=50,b=80,t=100,pad=5),paper_bgcolor='rgb(243, 243, 243)',
                         plot_bgcolor='rgb(243, 243, 243)',
                         annotations=annotation
                         )


        # labeling the bar net worht


    plotdiv = plot(fig, output_type = 'div', include_plotlyjs = False,show_link=False)
    return plotdiv
