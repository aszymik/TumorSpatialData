import pandas as pd
import os
import plotly.express as px
import streamlit as st

from helper import *
from patient_statistics import *

# funckja przechodzaca przez wszystkich pacjentów
# za dwa tyg (tylko prezentacja)

file_path = "if_data/"

patients = []
for filename in os.listdir(file_path + 'IF1'):
    patients.append(filename[:4])

st.set_page_config(layout="wide")
st.title('Tumor spatial data analysis')
with st.sidebar:
    st.title("Select parameters")
    patient = st.sidebar.selectbox('Select a patient', patients)
    panel = st.sidebar.selectbox('Select a panel', ['IF1', 'IF2', 'IF3'])

data = get_panel(panel, patient)

cell_types = ['all'] + list(data.celltype.unique())
cell_types_opt = st.multiselect(
    'Choose cell types',
    cell_types,
    default = 'all'
)

# if 'all' is selected, use all cell types; otherwise, use the selected ones
if 'all' in cell_types_opt:
    data_filtered = data
else:
    data_filtered = data[data.celltype.isin(cell_types_opt)]

fig = px.scatter(data_filtered, x='nucleus.x', y='nucleus.y', color="celltype", color_discrete_map=IF1_cell_mapping, title=f'Patient {patient} {panel} panel')
fig.update_layout(autosize=False, width=800, height=600)
st.plotly_chart(fig)

st.markdown(f'### Patient {patient} TLSs')

# Wizualizacja TLSów
TLS_graph = patient_TLSs_plot(data)
st.plotly_chart(TLS_graph)

# Bar plot TLSów pacjenta
bar_plot = patient_bar_plot(patient)
st.plotly_chart(bar_plot)

# Klastrowanie TLSów wszystkich pacjentów
st.markdown('### TLS differentiation')
all_patients_clusters = all_patients_clusters_plot()
st.plotly_chart(all_patients_clusters, use_container_width=True)


