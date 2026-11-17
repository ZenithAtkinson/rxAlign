import csv

pathlaptop=r"C:\Users\default.DESKTOP-UOBTTU9\Bioinformatics Lesgo\rxAlign\data\regimen_cancer_level_dataset.csv"
#pathhome=None

with open(pathlaptop, newline="") as f:
    drugs_drug = {} # dict init
    for row in csv.DictReader(f):
        
        if row["redcap_ca_index"] != "Yes": continue # skpping non-lungcancer stuff but can easily be removed if we want this in
        print(row["record_id"], row["regimen_number_within_cancer"], row["drugs_drug_1"]) #only printing first drug in each regimen

        # labeling. # lines 5+ grouped: few patients reach them (197 of 1,283), so separate columns would be mostly empty (adds noise? But we can go up to 8 or whatever)
        line = int(row["regimen_number_within_cancer"]) # turning it into ints for comparisons
        if line >= 5:
            label = "L5+"
        else:
            label = f"L{line}" 
        
        for n in range(1,6): # this is inefficient because not every one has 5
            drug = row[f"drugs_drug_{n}"]
            print("ze drug name is " + drug)
            if drug is "":
                continue
            drug = drug.split("(")[0] # removing the extra synonyms

            drugs_drug.setdefault(row["record_id"], set()).add(f"{drug}@{label}") # setdefault just maked an empty set for the first time, and .add just puts the given label in.

    print(len(drugs_drug))
    print(drugs_drug["GENIE-DFCI-000142"]) #returns order for this patient: {'Docetaxel@L5+', 'Investigational Drug@L4', 'Crizotinib@L5+', 'Cisplatin@L1', 'Pemetrexed Disodium@L3', 'Etoposide@L1', 'Investigational Drug@L2', 'Carboplatin@L3'}
        # the @ saves the order
    
    #making a new matrix using first matrix but adding this one:
    #Making all the labels unique so it matches across all patients
    all_labels = set()
    for labels in drugs_drug.values():
        all_labels.update(labels)
    all_labels = sorted(all_labels)

    # full 0/1 matrix
    patient_ids = list(drugs_drug.keys())
    new_matrix = []
    for pid in patient_ids:
        new_matrix.append([1 if label in drugs_drug[pid] else 0 for label in all_labels])


    # time to save it
    import os
    import pandas as pd

    df = pd.DataFrame(new_matrix, index=patient_ids, columns=all_labels)
    df.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "newer_timeline_matrix.csv"))
