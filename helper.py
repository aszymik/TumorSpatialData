import os
import pandas as pd

FILE_PATH = "if_data/"


def standardize_phenotype(phenotype):
    markers = set()
    start = 0
    for i in range(len(phenotype)):
        if phenotype[i] == '-' or phenotype[i] == '+':
            markers.add(phenotype[start:i + 1])
            start = i + 1
    return ''.join(sorted(markers))


def get_panel(panel, patient):
    mapping = pd.read_csv(f'dicts/{panel}_phen_to_cell_mapping.csv', sep=',', header=0)
    mapping['phenotype'] = mapping['phenotype'].apply(lambda x: standardize_phenotype(x))

    data = pd.read_csv(FILE_PATH + f"{panel}/" + patient + f"_{panel}.csv", sep=',', header=0, index_col=0)
    data['phenotype'] = data['phenotype'].apply(lambda x: standardize_phenotype(x))
    return data.merge(mapping, on='phenotype', how='left')


def get_all_patients(panel):

    dir_path = FILE_PATH + f"{panel}/"
    patients = []

    for filename in os.listdir(dir_path):
        if filename.endswith(f'_{panel}.csv'):
            patient = filename.split('_')[0]
            patients.append(patient)

    return patients


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

