
# graf tylko nad wybranym typem komórek (np. CK)
# chcemy sklastrowac wektory udzialu procentowego poszczegolnych typów komórek

# elementy na grafie TLS – szukamy sąsiadów w ogólnym grafie
# potem mozemy zobaczyc ile u naszych pacjentów jest danego klastra, jakies ciekawe podsumowanie


# kandydaci to teraz zbior krotek z id komórki
# wsrod kazdego kandydata (o jakiejs sensownej dlugosci) chcemy go zwizualizowac 
# i obliczyc procentowy udzial poszczegolnych typow komorek w tym klastrze
# moze z kazdego zrobic graf i zwizualizowac je na jednym obrazku?        


# Wizualizacja -- na razie trochę nie mam pomysłu
# x_values = []
# y_values = []
# colors = []

# for component in candidates:
#     for cell_id in component:
#         x, y = G_all.nodes[cell_id]['pos']
#         x_values.append(x)
#         y_values.append(y)
#         colors.append(IF1_cell_mapping[G_all.nodes[cell_id]['celltype']])

# fig = go.Figure(data=go.Scattergl(
#     x = x_values,
#     y = y_values,
#     mode = 'markers',
#     marker = dict(
#         color = colors,
#         size = 10
#     )
# ))

# fig.show()

# Wykresy dla kazdego pacjenta osobno
"""
fig, axes = plt.subplots(3, 6, figsize=(12, 6))
axes = axes.flatten()

for i, patient in enumerate(patients):
    ax = axes[i]
    
    patient_df = df[df['patient'] == patient]
    patient_df = patient_df.drop(columns=['patient', 'component_size'])
    patient_df.plot(kind='barh', stacked=True, ax=ax)

    if i == 0:
        handles, labels = ax.get_legend_handles_labels()

    ax.legend().remove()
    ax.set_title(f'Patient {patient}')
    adjust_axes(ax)

fig.legend(handles, labels, loc='upper center', ncol=9, fontsize=8, bbox_to_anchor=(0.5, 0.96))
fig.supxlabel('Cell percentage')
fig.supylabel('Clusters')
fig.suptitle('Cell distribution in TLS candidates', y=0.98)

plt.tight_layout()
plt.show()
"""