import pandas as pd
import os
import plotly.express as px
import streamlit as st

from helper import *
from patient_statistics import *

file_path = "if_data/"

patients = []
for filename in os.listdir(file_path + 'IF1'):
    patients.append(filename[:4])

st.set_page_config(layout="wide")
st.title('Tumor immune environment analysis')
with st.sidebar:
    st.title("Select parameters")
    patient = st.sidebar.selectbox('Select a patient', patients)
    panel = st.sidebar.selectbox('Select a panel', ['IF1', 'IF2', 'IF3'])

st.markdown('### Spatial distribution of cell types in chosen tissue fragment')

data = get_panel(panel, patient)
df = get_panel('IF1', patient) if panel != 'IF1' else data

cell_types = ['all'] + list(data.celltype.unique())
cell_types_opt = st.multiselect(
    'Choose cell types',
    cell_types,
    default = 'all'
)

if 'all' in cell_types_opt:  # if 'all' is selected, use all cell types; otherwise, use the selected ones
    data_filtered = data
else:
    data_filtered = data[data.celltype.isin(cell_types_opt)]

fig = px.scatter(data_filtered, x='nucleus.x', y='nucleus.y', opacity=1.0, color="celltype", color_discrete_map=IF1_cell_mapping, title=f'Patient {patient} {panel} panel visualisation')
fig.update_layout(autosize=False, width=800, height=600)
st.plotly_chart(fig)

# Procentowy udział komórek w tkance
cell_counts = data_filtered['celltype'].value_counts().reset_index()
cell_counts.columns = ['celltype', 'count']
tissue_bar = px.bar(cell_counts, x='celltype', y='count', color='celltype', color_discrete_map=IF1_cell_mapping, title='Cell type counts in the tissue')
tissue_bar.update_layout(autosize=False, width=800, height=600)
st.plotly_chart(tissue_bar)


st.markdown('### Tertiary Lymphoid Structures')
st.markdown('Tertiary Lymphoid Structures (TLS) are crucial components in the spatial arrangement of tumor tissues. They are a form of lymphoid tissue that emerge in non-lymphoid organs, such as tumors, due to persistent immune stimulation. TLS house various immune cells, including T cells, B cells and dendritic cells. They play a significant role in local immune responses and have been linked to better patient prognosis in certain types of cancer.')

# Sąsiedztwo kolejnych B-celli
st.markdown('#### B-cell neighborhood across tissue')
st.caption('Cells are sorted in ascending order by nucleus position.')
b_cell_neighbors = analyse_bcell_neighborhood(df)
st.plotly_chart(b_cell_neighbors)

# Wizualizacja TLSów
st.markdown(f'#### Patient {patient} TLS analysis')
TLS_graph = patient_TLS_plot(df)
st.plotly_chart(TLS_graph)

# Bar plot TLSów pacjenta
bar_plot = patient_bar_plot(patient)
st.plotly_chart(bar_plot)

# Klastrowanie TLSów wszystkich pacjentów
st.markdown('#### TLS differentiation across all patients')
all_patients_clusters = all_patients_clusters_plot()
st.plotly_chart(all_patients_clusters, use_container_width=True)


