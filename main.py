import networkx as nx
import plotly.graph_objects as go
from sklearn import neighbors
from scipy.spatial import KDTree
from helper import get_panel
from sklearn.metrics import pairwise_distances


patient = "1107"
df = get_panel('IF1', patient)
df.set_index('cell.ID')


def graph_by_cell_type(df, cell_types=None):
    """Tworzy graf sąsiedztwa dla wybranych typów komórek"""

    if cell_types is not None:
        df = df[df['celltype'].isin(cell_types)]

    adj_matrix = neighbors.radius_neighbors_graph(df[['nucleus.x', 'nucleus.y']], radius=30)  # macierz sąsiedztwa

    G = nx.Graph()  # graf z danymi o komórkach
    for i, (cell_id, cell_type, x, y) in enumerate(zip(df['cell.ID'], df['celltype'], df['nucleus.x'], df['nucleus.y'])):
        G.add_node(cell_id)  # wierzchołki to cell_id
        G.nodes[cell_id]['celltype'] = cell_type
        G.nodes[cell_id]['pos'] = x, y

    # Krawędzie na podstawie macierzy sąsiedztwa
    for i, row in enumerate(adj_matrix.toarray()):
        for j, connected in enumerate(row):
            if connected:
                G.add_edge(df.iloc[i]['cell.ID'], df.iloc[j]['cell.ID'])

    return G


def connected_component_from_candidate(tree, TLS_candidate, distance=30):

    def DFS(tree, start_point, visited=None, distance=distance):
        if visited is None:
            visited = set()

        # Znajdź wszystkie punkty w określonej odległości od punktu startowego
        indices = tree.query_ball_point(start_point, distance)

        for index in indices:
            if index not in visited:
                visited.add(index)
                DFS(tree, df.iloc[index][['nucleus.x', 'nucleus.y']], visited, distance)

        return visited

    # start_point = df[df['cell.ID'].isin(TLS_candidate.nodes)]['nucleus.x', 'nucleus.y'].iloc[0]
    start_cell_id = next(iter(TLS_candidate))
    start_coordinates = df.loc[start_cell_id, ['nucleus.x', 'nucleus.y']].values

    # Teraz 'coordinates' to współrzędne punktu o danym 'cell_id' w KDTree
    component = DFS(tree, start_coordinates)

    print("Spójna składowa zawierająca punkt startowy to:", df.iloc[list(component)])

    return component


# KDTree wszystkich komórek
# dist_matrix = pairwise_distances(df[['nucleus.x', 'nucleus.y']])
tree = KDTree(df[['nucleus.x', 'nucleus.y']])

# Spójne składowe, będące kandydatami na TLS
G_TLS_candidates = graph_by_cell_type(df, ['Tcell', 'Bcell', 'BnTcell'])
TLS_components = nx.connected_components(G_TLS_candidates)

for candidate in TLS_components:
    if len(candidate) > 5:
        connected_component_from_candidate(tree, candidate)



"""
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
"""
