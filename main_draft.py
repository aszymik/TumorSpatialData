import networkx as nx
import plotly.graph_objects as go
from sklearn import neighbors
from scipy.spatial import KDTree
from helper import get_panel

patient = "1107"
df = get_panel('IF1', patient)


# TODO: KDTree zamiast zwykłego grafu bo to mega długo trwa

# graf tylko nad wybranym typem komórek (np. CK)
# chcemy sklastrowac wektory udzialu procentowego poszczegolnych typów komórek

# elementy na grafie TLS – szukamy sąsaidów w ogólnym grafie
# potem mozemy zobaczyc ile u naszych pacjentów jest danego klastra, jakies ciekawe podsumowanie


def graph_by_cell_type(df, cell_types=None):
    """Tworzy graf sąsiedztwa dla wybranych typów komórek"""

    if cell_types is not None:
        df = df[df['celltype'].isin(cell_types)]

    adj_matrix = neighbors.radius_neighbors_graph(df[['nucleus.x', 'nucleus.y']], radius=30)  # macierz sąsiedztwa

    G = nx.Graph()  # Create an empty graph

    # Dodajemy dane o komórkach
    for i, (cell_id, cell_type, x, y) in enumerate(zip(df['cell.ID'], df['celltype'], df['nucleus.x'], df['nucleus.y'])):
        G.add_node(cell_id)  # Add node with cell_id
        G.nodes[cell_id]['celltype'] = cell_type
        G.nodes[cell_id]['pos'] = x, y

    # Add edges based on adjacency matrix
    for i, row in enumerate(adj_matrix.toarray()):
        for j, connected in enumerate(row):
            if connected:
                G.add_edge(df.iloc[i]['cell.ID'], df.iloc[j]['cell.ID'])

    return G

    G = nx.from_scipy_sparse_array(adj_matrix)  # graf na podstawie macierzy

    # Dodajemy dane o komórkach
    for i, (cell_id, cell_type, x, y) in enumerate(zip(df['cell.ID'], df['celltype'], df['nucleus.x'], df['nucleus.y'])):
        G.nodes[i]['cell_id'] = cell_id
        G.nodes[i]['celltype'] = cell_type
        G.nodes[i]['pos'] = x, y

    return G


G_all = graph_by_cell_type(df)
all_connected_components = nx.connected_components(G_all)
G_TLS_candidates = graph_by_cell_type(df, ['Tcell', 'Bcell', 'BnTcell'])
TLS_components = nx.connected_components(G_TLS_candidates)

candidates = {}

# # dla każdej spójnej składowej w kandydatach na TLS
# for cell in G_TLS_candidates:
#     for connected_component in all_connected_components:
#         # sprawdzamy, czy jakikolwiek z sąsiadów należy do spójnej składowej
#         if cell in connected_component:
#         # if any(cell in connected_component for cell in candidate):
#             G_TLS_candidates.add_nodes_from(connected_component)


for candidate in TLS_components:
    if len(candidate) < 5:
        continue
    start_cell_id = next(iter(candidate))
    for component in all_connected_components:
        if start_cell_id in component:
            candidates.add(component)

print(candidates)

# for i in G_TLS_candidates.nodes:
#     cell_id = G_TLS_candidates.nodes[i]['cell_id']
#     print(f'{cell_id=}')
#     for connected_component in all_connected_components:
#         print(f'{connected_component=}')
#         # sprawdzamy, czy jakikolwiek z sąsiadów należy do spójnej składowej
#         if cell_id in connected_component.nodes['cell_id']:
#         # if any(cell in connected_component for cell in candidate):
#             G_TLS_candidates.add_nodes_from(connected_component)   


# cell_id_to_node = {data['cell_id']: node for node, data in G_all.nodes(data=True)}
# print(f'{cell_id_to_node=}')

# # dla każdej spójnej składowej w kandydatach na TLS
# for cell_id in G_TLS_candidates:
#     for connected_component in all_connected_components:
#         # sprawdzamy, czy jakikolwiek z sąsiadów należy do spójnej składowej
#         if cell_id_to_node[cell_id] in connected_component:
#             G_TLS_candidates.add_nodes_from(connected_component)      


print(f'{G_TLS_candidates=}')
G = G_TLS_candidates

# cell_type_color_map = {
#     'Tcell': 'red',
#     'Bcell': 'blue',
#     'BnTcell': 'green',
# }

# node_colors = [cell_type_color_map[G.nodes[i]['celltype']] for i in G.nodes]
# node_colors = [G.nodes[i]['celltype'] for i in G.nodes]

celltypes = list(set(G.nodes[i]['celltype'] for i in G.nodes))
celltype_to_number = {celltype: number for number, celltype in enumerate(celltypes)}

# Teraz możemy użyć tego słownika do kolorowania wierzchołków
node_colors = [celltype_to_number[G.nodes[i]['celltype']] for i in G.nodes]


edge_x = []
edge_y = []
for edge in G.edges():
    x0, y0 = G.nodes[edge[0]]['pos']
    x1, y1 = G.nodes[edge[1]]['pos']
    edge_x.append(x0)
    edge_x.append(x1)
    edge_x.append(None)
    edge_y.append(y0)
    edge_y.append(y1)
    edge_y.append(None)

edge_trace = go.Scatter(
    x=edge_x, y=edge_y,
    line=dict(width=0.5, color='#888'),
    hoverinfo='none',
    mode='lines')

node_x = []
node_y = []
for node in G.nodes():
    x, y = G.nodes[node]['pos']
    node_x.append(x)
    node_y.append(y)

node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode='markers',
    hoverinfo='text',
    marker=dict(
        showscale=True,
        colorscale='YlGnBu',
        reversescale=True,
        color=node_colors,
        size=10,
        colorbar=dict(
            thickness=15,
            title='Node Connections',
            xanchor='left',
            titleside='right'
        ),
        line_width=2))

# node_adjacencies = []
# node_text = []
# for node, adjacencies in enumerate(G.adjacency()):
#     node_adjacencies.append(len(adjacencies[1]))
#     node_text.append('# of connections: '+str(len(adjacencies[1])))

# node_trace.marker.color = node_adjacencies
# node_trace.text = node_text

fig = go.Figure(data=[edge_trace, node_trace],
             layout=go.Layout(
                title='<br>Network graph made with Python',
                titlefont_size=16,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )
fig.show()
