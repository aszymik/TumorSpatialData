import pandas as pd
import plotly.figure_factory as ff
import plotly.graph_objects as go
import plotly.express as px

from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import ward
from plotly.graph_objs import Figure
from helper import *
from main import *


def all_patients_clusters_plot() -> Figure:
    """Procentowy udział typów komórek w TLS-ach wszystkich pacjentów"""

    df = pd.read_csv('if_data/cell_type_percentages_in_TLS.tsv', sep='\t', index_col=0)
    df = df.drop('component_size', axis=1)

    # Dendrogram
    colors = [IF1_cell_mapping[celltype] for celltype in ('Macrophage', 'Tcell', 'other', 'Bcell', 'CD15+Tumor', 'CD15-Tumor')]
    fig = ff.create_dendrogram(df.drop('patient', axis=1), orientation='right', linkagefun=ward, labels=df.index, color_threshold=3.5, colorscale=colors)
    dendro_side = go.Figure(data=fig['data'])

    # Bar plot
    df = df.reindex(fig['layout']['yaxis']['ticktext'])
    df.index = df.index.astype(str)

    bar_side = go.Figure()
    for col in df.drop('patient', axis=1).columns:
        bar_side.add_trace(go.Bar(y=fig['layout']['yaxis']['tickvals'], x=df[col], name=col, orientation='h',
                                  marker=dict(color=IF1_cell_mapping[col])))

    # Łączymy dendrogram i bar plot
    fig = make_subplots(rows=1, cols=2, shared_yaxes=True, horizontal_spacing=0)

    # Dodajemy dane
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
    fig.layout.title = 'Cell type distribution in TLS candidates'
    fig.layout.plot_bgcolor = 'rgba(0,0,0,0)'

    fig.update_layout(autosize=False, width=900, height=600)
    return fig


def patient_bar_plot(patient: str) -> Figure:
    """Procentowy udział typów komórek w TLS-ach danego pacjenta"""

    df = pd.read_csv('if_data/cell_type_percentages_in_TLS.tsv', sep='\t', index_col=0)
    df = df[df['patient'] == int(patient)]
    df = df.drop('component_size', axis=1)

    bar_plot = go.Figure()
    for col in df.drop('patient', axis=1).columns:
        bar_plot.add_trace(go.Bar(y=df.index, x=df[col], name=col, orientation='h',
                              marker=dict(color=IF1_cell_mapping[col])))

    bar_plot.layout.xaxis.title = 'Cell percentage'
    bar_plot.layout.xaxis.range = [0, 1]
    bar_plot.layout.barmode = 'stack'
    bar_plot.layout.title = f'Cell type distribution in patient {patient} TLS'

    bar_plot.update_layout(autosize=False, width=900, height=600)
    return bar_plot


def patient_TLS_plot(df: pd.DataFrame) -> Figure:
    """Wizualizuje TLS-y pacjenta na tkance"""

    _, candidates = TLS_candidates(df)
    candidate_nodes = set()
    for candidate in candidates:
        candidate_nodes.update(candidate)  

    df = df.copy()
    df.loc[~df['cell.ID'].isin(candidate_nodes), 'celltype'] = 'not a TLS'  # kolorujemy wg typu komórki tylko elementy z TLSów
    fig = px.scatter(df, x='nucleus.x', y='nucleus.y', opacity=0.5, color='celltype', color_discrete_map=IF1_cell_mapping, title='Spatial distribution of TLS in the tissue')
    fig.update_traces(marker_size=3)
    fig.update_layout(autosize=False, width=800, height=600)
    return fig


def analyse_bcell_neighborhood(df: pd.DataFrame) -> Figure:
    """Analizuje typy komórek w sąsiedztwie B-celli i wizualizuje je"""

    G_all = graph_by_cell_type(df)
    b_cell_df = df[df['celltype'] == 'Bcell']

    neighbors_data = []
    for _, row in b_cell_df.iterrows():
        # Zliczamy kazdy typ komórki
        celltype_counts = {cell_type: 0 for cell_type in CELL_TYPES}
        cell_neighbors = list(G_all.neighbors(row['cell.ID']))
        for neighbor in cell_neighbors:
            celltype = G_all.nodes[neighbor]['celltype']
            celltype_counts[celltype] += 1

        # Obliczamy procentowy udział
        num_neighbors = len(cell_neighbors)
        neighbors_dict = {cell_type: (count / num_neighbors) for cell_type, count in celltype_counts.items()}
        
        # Dodajemy 'nucleus.x' i 'nucleus.y'  z oryginalnej ramki danych
        neighbors_dict['cell.ID'] = row['cell.ID']
        neighbors_dict['nucleus.x'] = row['nucleus.x']
        neighbors_dict['nucleus.y'] = row['nucleus.y']

        neighbors_data.append(neighbors_dict)

    # Tworzymy DataFrame
    neighbors_df = pd.DataFrame(neighbors_data, columns=['cell.ID', 'nucleus.x', 'nucleus.y', 'other', 'CD15-Tumor', 'CD15+Tumor', 'Tcell', 'Bcell', 'BnTcell', 'Macrophage', 'Neutrophil', 'DC'])
    neighbors_df = neighbors_df.sort_values(by=['nucleus.x', 'nucleus.y']).reset_index()
    neighbors_df = neighbors_df.drop(columns=['index'])


    bar_plot = go.Figure()
    for col in neighbors_df.drop(columns=['cell.ID', 'nucleus.x', 'nucleus.y']).columns:
        bar_plot.add_trace(go.Bar(y=neighbors_df.index, x=neighbors_df[col], name=col, orientation='h',
                              marker=dict(color=IF1_cell_mapping[col])))

    bar_plot.layout.xaxis.title = 'Cell percentage'
    bar_plot.layout.xaxis.range = [0, 1]
    bar_plot.layout.barmode = 'stack'
    bar_plot.layout.title = f'Cell type distribution across B-cell neighbors'
    bar_plot.update_layout(autosize=False, width=700, height=500)

    return bar_plot
