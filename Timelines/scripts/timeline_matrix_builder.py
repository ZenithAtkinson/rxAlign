import os
import sys

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

all_drugs = {}
for patient in patient_drugs:
    for drug in patient_drugs[patient]:
        for d in drugs:
            print(d)
