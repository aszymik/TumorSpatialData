import networkx as nx
import plotly.graph_objects as go
from sklearn import neighbors
from helper import get_panel

patient = "1107"
df = get_panel('IF1', patient)
celltypes = df['celltype'].unique()

def rgb_to_hex(rgb):
    return '#' + '%02x%02x%02x' % rgb

IF1_cell_mapping = {"other": rgb_to_hex((190, 190, 190)), 
                    "CD15+Tumor": rgb_to_hex((73, 176, 248)),
                    "CD15-Tumor": rgb_to_hex((138, 79, 45)),
                    "Tcell": rgb_to_hex((235, 74, 148)),
                    "Bcell": rgb_to_hex((204, 49, 31)),
                    "BnTcell": rgb_to_hex((236, 95, 42)),
                    "Neutrophil": rgb_to_hex((0, 40, 245)),
                    "Macrophage": rgb_to_hex((97, 209, 62)),
                    "DC": rgb_to_hex((49, 113, 30))}


# graf tylko nad wybranym typem komórek (np. CK)
# chcemy sklastrowac wektory udzialu procentowego poszczegolnych typów komórek

# elementy na grafie TLS – szukamy sąsiadów w ogólnym grafie
# potem mozemy zobaczyc ile u naszych pacjentów jest danego klastra, jakies ciekawe podsumowanie


def graph_by_cell_type(df, cell_types=None):
    """Tworzy graf sąsiedztwa dla wybranych typów komórek"""

    if cell_types is not None:
        df = df[df['celltype'].isin(cell_types)]

    adj_matrix = neighbors.radius_neighbors_graph(df[['nucleus.x', 'nucleus.y']], radius=30)  # macierz sąsiedztwa  
    G = nx.from_scipy_sparse_array(adj_matrix)

    index_to_cell_id = df['cell.ID'].to_dict()
    G = nx.relabel_nodes(G, index_to_cell_id)

    # Dodajemy dane o komórkach
    for cell_id, cell_type, x, y in zip(df['cell.ID'], df['celltype'], df['nucleus.x'], df['nucleus.y']):
        if cell_id in G.nodes:
            G.nodes[cell_id]['celltype'] = cell_type
            G.nodes[cell_id]['pos'] = x, y

    return G


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
        if start_cell_id in component:
            candidates.add(tuple(component))

# kandydaci to teraz zbior krotek z id komórki
# wsrod kazdego kandydata (o jakiejs sensownej dlugosci) chcemy go zwizualizowac 
# i obliczyc procentowy udzial poszczegolnych typow komorek w tym klastrze
# moze z kazdego zrobic graf i zwizualizowac je na jednym obrazku?        



