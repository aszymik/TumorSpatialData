import pandas as pd
import plotly.figure_factory as ff
import plotly.graph_objects as go
import plotly.express as px

from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import ward
from helper import *
from main import TLS_candidates

# Przykładowa lista kolorów


def all_patients_clusters_plot():
    df = pd.read_csv('if_data/cell_type_percentages_in_TLSs.tsv', sep='\t', index_col=0)
    df = df.drop('component_size', axis=1)

    # Dendrogram
    colors = [IF1_cell_mapping[celltype] for celltype in ('Macrophage', 'CD15-Tumor', 'other', 'Bcell', 'CD15+Tumor', 'Tcell')]
    fig = ff.create_dendrogram(df.drop('patient', axis=1), orientation='right', linkagefun=ward, labels=df.index, color_threshold=1.5, colorscale=colors)
    dendro_side = go.Figure(data=fig['data'])

    # Bar plot
    df = df.reindex(fig['layout']['yaxis']['ticktext'])
    df.index = df.index.astype(str)

    bar_side = go.Figure()
    for col in df.drop('patient', axis=1).columns:
        bar_side.add_trace(go.Bar(y=fig['layout']['yaxis']['tickvals'], x=df[col], name=col, orientation='h',
                                  marker=dict(color=IF1_cell_mapping[col])))

    # Combine dendrogram and bar plot
    fig = make_subplots(rows=1, cols=2, shared_yaxes=True, horizontal_spacing=0)

    # Add data
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i].showlegend = False
        fig.add_trace(dendro_side['data'][i], row=1, col=1)

    for i in range(len(bar_side['data'])):
        fig.add_trace(bar_side['data'][i], row=1, col=2)

    fig.layout.xaxis.showticklabels = False
    fig.layout.yaxis.showticklabels = False
    fig.layout.xaxis2.title = 'Cell percentage'
    fig.layout.xaxis2.range = [0, 1]
    fig.layout.barmode = 'stack'
    fig.layout.title = 'Cell distribution in TLS candidates'
    fig.layout.plot_bgcolor = 'rgba(0,0,0,0)'

    fig.update_layout(autosize=False, width=900, height=600)
    return fig


def patient_bar_plot(patient):
    df = pd.read_csv('if_data/cell_type_percentages_in_TLSs.tsv', sep='\t', index_col=0)
    df = df[df['patient'] == int(patient)]
    df = df.drop('component_size', axis=1)

    bar_plot = go.Figure()
    for col in df.drop('patient', axis=1).columns:
        bar_plot.add_trace(go.Bar(y=df.index, x=df[col], name=col, orientation='h',
                              marker=dict(color=IF1_cell_mapping[col])))

    bar_plot.layout.xaxis.title = 'Cell percentage'
    bar_plot.layout.xaxis.range = [0, 1]
    bar_plot.layout.barmode = 'stack'
    bar_plot.layout.title = f'Cell distribution in patient {patient} TLSs'
    bar_plot.layout.plot_bgcolor = 'rgba(0,0,0,0)'

    bar_plot.update_layout(autosize=False, width=900, height=600)
    return bar_plot


def patient_TLSs_plot(df):
    G_all, candidates = TLS_candidates(df)

    node_x = []
    node_y = []
    node_labels = []
    for node in G_all.nodes():
        x, y = G_all.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)
        if any(node in candidate for candidate in candidates):
            node_labels.append(G_all.nodes[node]['celltype'])  # koloruj według typu komórki
        else:
            node_labels.append('not applicable')  # koloruj szaro pozostałe wierzchołki


    fig = px.scatter(x=node_x, y=node_y, color=node_labels, color_discrete_map=IF1_cell_mapping, title=f'')
    fig.update_layout(autosize=False, width=800, height=600)
    return fig