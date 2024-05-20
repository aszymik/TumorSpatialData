import networkx as nx
from sklearn import neighbors
from typing import List, Iterable, Tuple
from helper import *

CELL_TYPES = ['other', 'CD15-Tumor', 'CD15+Tumor', 'Tcell', 'Bcell', 'BnTcell', 'Macrophage', 'Neutrophil', 'DC']
PATIENTS = get_all_patients('IF1')


def graph_by_cell_type(df: pd.DataFrame, cell_types: Iterable[str] = None, radius: int = 30) -> nx.Graph:
    """Tworzy graf sąsiedztwa dla wybranych typów komórek."""

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


def TLS_candidates(df: pd.DataFrame, component_min: int = 20) -> Tuple[nx.Graph, List[set]]:
    """Znajduje spójne składowe oparte na B i T-cellach i dodaje do nich sąsiadów komórek z pełnego grafu.
       Zwraca pełny graf sąsiedztwa oraz listę zbiorów – kandydatów na TLS."""

    G_all = graph_by_cell_type(df)
    G_TLS_candidates = graph_by_cell_type(df, ['Tcell', 'Bcell', 'BnTcell'])
    TLS_components = list(nx.connected_components(G_TLS_candidates))
    TLS_components = list(filter(lambda x: len(x) >= component_min, TLS_components))  # filtrujemy po długości

    for i in range(len(TLS_components)):
        component = TLS_components[i]
        neighbours_set = set()
        for cell in component:
            neighbors = list(G_all.neighbors(cell))
            neighbours_set.update(neighbors)  
        component.update(neighbours_set)  # dodajemy sąsiadow z ogólnego grafu do składowej

    return G_all, TLS_components


def cell_types_in_patient_TLSs(patient: str, G_all: nx.Graph, candidates: List[set]) -> pd.DataFrame:
    """Tworzy DataFrame z procentowym udziałem poszczególnych typów komórek w kandydatach na TLS u danego pacjenta"""

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

    return pd.DataFrame(data)


def main() -> None:
    """Zbiera dane TLS-ów wszystkich pacjentów i zapisuje je do pliku, aby mieć łatwy dostęp do danych"""

    df = pd.DataFrame(columns=['patient', 'component_size', 'other', 'CD15-Tumor', 'CD15+Tumor', 'Tcell', 'Bcell', 'BnTcell', 'Macrophage', 'Neutrophil', 'DC'])

    for patient in PATIENTS:
        patient_df = get_panel('IF1', patient)

        # Graf sąsiedztwa i kandydaci na TLS (spójne składowe)
        G_all, candidates = TLS_candidates(patient_df)  

        # Procentowy udział komórek w TLS
        patient_cell_percentage = cell_types_in_patient_TLSs(patient, G_all, candidates)  

        # Dodajemy do ogólnej ramki danych
        df = pd.concat((df, patient_cell_percentage), ignore_index=True)

    df.to_csv('if_data/cell_type_percentages_in_TLS.tsv', sep='\t')


if __name__ == '__main__':
    main()
