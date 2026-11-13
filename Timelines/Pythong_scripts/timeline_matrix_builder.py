import os
import sys
import pandas as pd
import numpy as np

CLASSES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Classes")
sys.path.insert(0, CLASSES_DIR)

from load_timelines import load_timelines

pathlaptop=r"C:\Users\default.DESKTOP-UOBTTU9\Bioinformatics Lesgo\rxAlign\data\regimen_cancer_level_dataset.csv"
#pathhome=None
ze_timelines = load_timelines(pathlaptop)

patient_drugs = {}

for patient_id, timeline in ze_timelines.items():
    #print(ze_timelines[i].event_list[0].treatments)
    drugs = set()
    for event in timeline.event_list:
        drugs.update(event.treatments.treatment)
    patient_drugs[patient_id] = drugs

#print(len(patient_drugs))
#first = list(patient_drugs.keys())[0]
#print(first)

# put every drug into a set so we have a unique list
all_drugs = set()
for drugs in patient_drugs.values():
    all_drugs.update(drugs)
all_drugs = sorted(all_drugs)
#print(all_drugs)

patient_ids = list(patient_drugs.keys())

# build the matrix (where rows are ID's and the columns are boolean for each given drug). rows are ID's because thats what sctlearn is looking for
matrix = []
for drug in all_drugs:
    row = []
    for pid in patient_ids:
        row.append(1 if drug in patient_drugs[pid] else 0)
    matrix.append(row)

#turning it into a dataframe so its easier to work with
df = pd.DataFrame(matrix, index=all_drugs, columns=patient_ids)
df = df.T

df.to_csv("drug_matrix.csv")

