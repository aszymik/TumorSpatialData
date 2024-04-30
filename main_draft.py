import networkx as nx
import plotly.graph_objects as go
from sklearn import neighbors
from helper import get_panel

patient = "1107"
df = get_panel('IF1', patient)
celltypes = df['celltype'].unique()


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


print(f'{all_connected_components=}')
G = G_TLS_candidates


# node_colors = [cell_type_color_map[G.nodes[i]['celltype']] for i in G.nodes]
# node_colors = [G.nodes[i]['celltype'] for i in G.nodes]

# celltypes = list(set(G.nodes[i]['celltype'] for i in G.nodes))
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
