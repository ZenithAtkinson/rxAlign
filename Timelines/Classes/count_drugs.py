# Prints every unique value in the drugs_drug_1 ... drugs_drug_5 columns and how many times each appears.
# Values are counted exactly as they appear in the file (no cleaning or combining).
# Usage: python count_drugs.py [path_to_csv]
import csv
import os
import sys
from collections import Counter

default_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "regimen_cancer_level_dataset.csv")
path = sys.argv[1] if len(sys.argv) > 1 else default_path

counts = Counter()
with open(path, newline="") as f:
    for row in csv.DictReader(f):
        for n in range(1, 6):
            drug = row[f"drugs_drug_{n}"]
            if drug == "":      # this regimen has fewer than n drugs
                continue
            counts[drug] += 1

for drug, count in counts.most_common():
    print(f"{count}\t{drug}")

print(f"\n{len(counts)} unique drug names, {sum(counts.values())} total")
