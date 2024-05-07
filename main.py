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

    index_to_cell_id = df['cell.ID'].to_dict()
    G = nx.relabel_nodes(G, index_to_cell_id)

    # Dodajemy dane o komórkach
    for cell_id, cell_type in zip(df['cell.ID'], df['celltype']):
        if cell_id in G.nodes:
            G.nodes[cell_id]['celltype'] = cell_type
        else:   
            G.add_node(cell_id, celltype=cell_type)
    return G


def TLS_candidates(df, component_min=20):
    """Zwraca spójne składowe grafu wzbogacone w komórki odpornościowe"""
    G_all = graph_by_cell_type(df)
    all_connected_components = list(nx.connected_components(G_all))
    G_TLS_candidates = graph_by_cell_type(df, ['Tcell', 'Bcell', 'BnTcell'])
    TLS_components = list(nx.connected_components(G_TLS_candidates))

    candidates = set()
    for candidate in TLS_components:
        if len(candidate) < component_min:
            continue
        for component in all_connected_components:
            if len(component) >= component_min:
                if any(cell_id in component for cell_id in candidate):
                    candidates.add(tuple(component))

    # print(f'{len(candidates)=}')
    return G_all, candidates


def cell_types_in_patient_TLSs(patient, G_all, candidates):
    """Zwraca DataFrame z procentowym udziałem poszczególnych typów komórek w kandydatach na TLS u danego pacjenta"""

    data = []
    for component in candidates:
        # Zliczamy kazdy typ komórek
        celltype_counts = {cell_type: 0 for cell_type in CELL_TYPES}
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
    df = pd.DataFrame(data, columns=['patient', 'component_size', 'other', 'CD15-Tumor', 'CD15+Tumor', 'Tcell', 'Bcell', 'BnTcell', 'Macrophage', 'Neutrophil', 'DC'])
    return df


if __name__ == '__main__':
    PATIENTS = get_all_patients('IF1')
    df = pd.DataFrame(columns=['patient', 'component_size', 'other', 'CD15-Tumor', 'CD15+Tumor', 'Tcell', 'Bcell', 'BnTcell', 'Macrophage', 'Neutrophil', 'DC'])

    for patient in PATIENTS:
        print(f'PATIENT {patient}')
        patient_df = get_panel('IF1', patient)  # wszystkie dane pacjenta
        G_all, candidates = TLS_candidates(patient_df)  # graf sąsiedztwa i kandydaci na TLS (spójne składowe)
        patient_cell_percentage = cell_types_in_patient_TLSs(patient, G_all, candidates)  # procentowy udział komórek w TLS
        df = pd.concat((df, patient_cell_percentage), ignore_index=True)  # dodajemy do ogólnej ramki danych

    df.to_csv('if_data/cell_type_percentages_in_TLSs.tsv', sep='\t')

