import networkx as nx
import plotly.graph_objects as go
from sklearn import neighbors
from helper import *

# TODO: zmień kolejność
cell_types = ['other', 'CD15-Tumor', 'Tcell', 'Bcell', 'Macrophage', 'BnTcell', 'Neutrophil', 'CD15+Tumor', 'DC']


def graph_by_cell_type(df, cell_types=None):
    """Tworzy graf sąsiedztwa dla wybranych typów komórek"""

    if cell_types is not None:
        df = df[df['celltype'].isin(cell_types)]

    adj_matrix = neighbors.radius_neighbors_graph(df[['nucleus.x', 'nucleus.y']], radius=30, include_self=True)  # macierz sąsiedztwa  
    G = nx.from_scipy_sparse_array(adj_matrix)

    index_to_cell_id = df['cell.ID'].to_dict()
    G = nx.relabel_nodes(G, index_to_cell_id)

    # Dodajemy dane o komórkach
    for cell_id, cell_type, x, y in zip(df['cell.ID'], df['celltype'], df['nucleus.x'], df['nucleus.y']):
        if cell_id in G.nodes:
            G.nodes[cell_id]['celltype'] = cell_type
            G.nodes[cell_id]['pos'] = x, y

    return G


def TLS_candidates(df):
    """Zwraca spójne składowe grafu wzbogacone w komórki odpornościowe"""
    G_all = graph_by_cell_type(df)
    all_connected_components = list(nx.connected_components(G_all))
    G_TLS_candidates = graph_by_cell_type(df, ['Tcell', 'Bcell', 'BnTcell'])
    TLS_components = nx.connected_components(G_TLS_candidates)

    candidates = set()
    for candidate in TLS_components:
        if len(candidate) < 5:
            continue
        start_cell_id = next(iter(candidate))
        for component in all_connected_components:
            if len(component) > 20:
                if start_cell_id in component:
                    candidates.add(tuple(component))

    return G_all, candidates



# def cell_types_in_patient_TLSs(patient):
#     """Zwraca wektory procentowego udziału poszczególnych typów komórek w kandydatach na TLS u danego pacjenta"""

#     # TODO: przemyśleć czy nie lepszy będzie dataframe albo numpy array

#     df = get_panel('IF1', patient)
#     G_all, candidates = TLS_candidates(df)

#     percentage_list = []
#     for component in candidates:
#         # Zliczamy kazdy typ komórek
#         celltype_counts = {cell_type: 0 for cell_type in cell_types}
        
#         for cell_id in component:
#             celltype = G_all.nodes[cell_id]['celltype']
#             if celltype in celltype_counts:
#                 celltype_counts[celltype] += 1
        
#         # Udział każdego typu komórki
#         total_cells = len(component)
#         celltype_percentage = {cell_type: (count / total_cells) for cell_type, count in celltype_counts.items()}
#         celltype_percentage['component_size'] = total_cells

#         percentage_list.append(celltype_percentage)

#     return percentage_list    


def cell_types_in_patient_TLSs(patient):
    """Zwraca DataFrame z procentowym udziałem poszczególnych typów komórek w kandydatach na TLS u danego pacjenta"""

    df = get_panel('IF1', patient)
    G_all, candidates = TLS_candidates(df)

    data = []
    for component in candidates:
        # Zliczamy kazdy typ komórek
        celltype_counts = {cell_type: 0 for cell_type in cell_types}
        for cell_id in component:
            celltype = G_all.nodes[cell_id]['celltype']
            if celltype in celltype_counts:
                celltype_counts[celltype] += 1
        
        # Udział każdego typu komórki
        total_cells = len(component)
        celltype_percentage = {cell_type: (count / total_cells) for cell_type, count in celltype_counts.items()}
        celltype_percentage['component_size'] = total_cells
        celltype_percentage['patient'] = patient

        data.append(celltype_percentage)

    # Tworzymy DataFrame
    df = pd.DataFrame(data, columns=['patient', 'component_size', 'other', 'CD15-Tumor', 'Tcell', 'Bcell', 'Macrophage', 'BnTcell', 'Neutrophil', 'CD15+Tumor', 'DC'])
    return df


patients = get_all_patients('IF1')
df = pd.DataFrame(columns=['patient', 'component_size', 'other', 'CD15-Tumor', 'Tcell', 'Bcell', 'Macrophage', 'BnTcell', 'Neutrophil', 'CD15+Tumor', 'DC'])


for patient in patients:
    print(f'PATIENT {patient}')
    patient_df = cell_types_in_patient_TLSs(patient)
    df = pd.concat((df, patient_df), ignore_index=True)
    print(df)

df.to_csv('if_data/cell_type_percentages_in_TLSs.tsv', sep='\t')

# klastrujemy niezaleznie od pacjenta, potem patrzymy ile u pacjenta jest kazdego z klastrów
