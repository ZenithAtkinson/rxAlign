import os
import sys

CLASSES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Classes")
sys.path.insert(0, CLASSES_DIR)

from load_timelines import load_timelines
pathlaptop=r"C:\Users\default.DESKTOP-UOBTTU9\Bioinformatics Lesgo\rxAlign\data\regimen_cancer_level_dataset.csv"
#pathhome=None
ze_timelines = load_timelines(r"C:\Users\default.DESKTOP-UOBTTU9\Bioinformatics Lesgo\rxAlign\data\regimen_cancer_level_dataset.csv")

for i in ze_timelines:
    print(ze_timelines[i].event_list[0].treatments)
