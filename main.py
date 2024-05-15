import networkx as nx
from sklearn import neighbors
from helper import *

CELL_TYPES = ['other', 'CD15-Tumor', 'CD15+Tumor', 'Tcell', 'Bcell', 'BnTcell', 'Macrophage', 'Neutrophil', 'DC']


def graph_by_cell_type(df, cell_types=None, radius=30):
    """Tworzy graf sąsiedztwa dla wybranych typów komórek"""

    if cell_types is not None:
        df = df[df['celltype'].isin(cell_types)]

    adj_matrix = neighbors.radius_neighbors_graph(df[['nucleus.x', 'nucleus.y']].values, radius=radius, include_self=True)  # macierz sąsiedztwa 
    G = nx.from_scipy_sparse_array(adj_matrix)
    index_to_cell_id = {i: id for i, id in enumerate(df['cell.ID'])}  # zmieniamy indeksy wierzchołków na id komórek
    G = nx.relabel_nodes(G, index_to_cell_id)

    # Dodajemy dane o komórkach
    for id, cell_type in zip(df['cell.ID'], df['celltype']):
        if id in G.nodes:
            G.nodes[id]['celltype'] = cell_type
        else:   
            G.add_node(id, celltype=cell_type)
    return G


def TLS_candidates(df, component_min=20):
    """Zwraca spójne składowe grafu wzbogacone w komórki odpornościowe"""
    G_all = graph_by_cell_type(df)
    all_connected_components = list(nx.connected_components(G_all))
    all_connected_components = list(filter(lambda x: len(x) >= component_min, all_connected_components))  # filtrujemy na podstawie długości
    # all_connected_components = [comp for comp in all_connected_comps if len(comp) >= component_min]    

    G_TLS_candidates = graph_by_cell_type(df, ['Tcell', 'Bcell', 'BnTcell'])
    TLS_components = list(nx.connected_components(G_TLS_candidates))

    TLS_components = list(filter(lambda x: len(x) >= component_min, TLS_components))
    # TLS_components = [comp for comp in TLS_comps if len(comp) >= component_min]   

    candidates = set()
    for candidate in TLS_components:
        for component in all_connected_components:
                if any(cell_id in component for cell_id in candidate):  # dodajemy sąsiadów z ogólnego grafu
                    candidates.add(tuple(component))
       
    return G_all, candidates

def cell_types_in_patient_TLSs(patient, G_all, candidates):
    """Zwraca DataFrame z procentowym udziałem poszczególnych typów komórek w kandydatach na TLS u danego pacjenta"""

    data = []
    for component in candidates:
        # Zliczamy kazdy typ komórek
        celltype_counts = {cell_type: 0 for cell_type in CELL_TYPES}
        for cell_id in component:
            celltype = G_all.nodes[cell_id]['celltype']
            celltype_counts[celltype] += 1

        # Udział każdego typu komórki
        total_cells = len(component)
        celltype_percentage = {cell_type: (count / total_cells) for cell_type, count in celltype_counts.items()}
        celltype_percentage['component_size'] = total_cells
        celltype_percentage['patient'] = patient
        data.append(celltype_percentage)

    return pd.DataFrame(data)  # zapisujemy do DataFrame


if __name__ == '__main__':
    PATIENTS = get_all_patients('IF1')
    df = pd.DataFrame(columns=['patient', 'component_size', 'other', 'CD15-Tumor', 'CD15+Tumor', 'Tcell', 'Bcell', 'BnTcell', 'Macrophage', 'Neutrophil', 'DC'])

    for patient in PATIENTS:
        print(f'PATIENT {patient}')
        patient_df = get_panel('IF1', patient)  # wszystkie dane pacjenta
        G_all, candidates = TLS_candidates(patient_df)  # graf sąsiedztwa i kandydaci na TLS (spójne składowe)
        patient_cell_percentage = cell_types_in_patient_TLSs(patient, G_all, candidates)  # procentowy udział komórek w TLS
        df = pd.concat((df, patient_cell_percentage), ignore_index=True)  # dodajemy do ogólnej ramki danych
        print(df)

    df.to_csv('if_data/cell_type_percentages_in_TLSs.tsv', sep='\t')
