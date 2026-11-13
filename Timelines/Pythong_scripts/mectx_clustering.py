# Python version of the MEC-TX method (R package "mectx", Dhrubo & Spakowicz 2026):
#   1. map each drug to a treatment class (Chemo, IO, Targeted, Hormone, Ancillary, Others)
#   2. put each patient on a 4-week time grid, starting at their first treatment, up to 5 years
#   3. binary feature for every (time bin, class): was the patient on that class in that bin?
#   4. drop near-constant features, PCA (no scaling, up to 50 PCs)
#   5. k-means on the PCA scores for every k in 3..20
# Differences from the R package: no Radiation/Small_Molecule classes (not in this data), and
# PCA is used for the 2D plot instead of UMAP.
# Run from Timelines/scripts. Outputs go to the current folder.
import csv
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score

pathlaptop = r"C:\Users\default.DESKTOP-UOBTTU9\Bioinformatics Lesgo\rxAlign"
#pathhome = r"C:\Users\Zenith\PhDCoding\rxAlign"
ROOT = pathlaptop
DATA_PATH = ROOT + r"\data\regimen_cancer_level_dataset.csv"
MED_CLASS_PATH = ROOT + r"\Code\data\med_class.tsv"

GRID_DAYS = 28              # 4-week bins (mectx grid_weeks = 4)
HORIZON_DAYS = 5 * 365.25   # 5 years (mectx horizon_years = 5)
MIN_FEATURE_VARIANCE = 0.01
N_PCS = 50
K_RANGE = range(3, 21)
FINAL_K = 6                 # same k as the Level 1 k-means, so the two can be compared
INDEX_CANCER_ONLY = True    # skip regimens for a patient's other cancers (redcap_ca_index == "No")

# ---------- 1. drug -> class ----------
# start from the previous students' med_class.tsv (Targeted / Chemo / ICB / Other)
CLASS_RENAME = {"Targeted": "Targeted", "Chemo": "Chemo", "ICB": "IO", "Other": "Others"}
drug_class = {}
with open(MED_CLASS_PATH, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        drug_class[row["Medication"].strip().upper()] = CLASS_RENAME[row["Class"]]

# drugs missing from med_class.tsv, plus hormone therapies (mectx has a separate Hormone class)
CLASS_OVERRIDES = {
    "Targeted": ["AFATINIB DIMALEATE", "CABOZANTINIB SMALATE", "DARATUMUMAB", "ENSARTINIB", "ENTRECTINIB",
                 "ERLOTINIB HYDROCHLORIDE", "ICOTINIB", "NAVITOCLAX", "NECITUMUMAB", "RITUXIMAB AND HYALURONIDASE HUMAN",
                 "ROCILETINIB", "RUXOLITINIB PHOSPHATE", "TRASTUZUMAB EMTANSINE"],
    "Chemo": ["ARSENIC TRIOXIDE", "BLEOMYCIN", "BUSULFAN", "CHLORAMBUCIL", "ETOPOSIDE PHOSPHATE", "LURBINECTEDIN",
              "MECHLORETHAMINE HYDROCHLORIDE", "PACLITAXEL LOADED POLYMERIC MICELLE", "PACLITAXEL POLIGLUMEX",
              "PACLITAXEL TREVATIDE", "PEGYLATED LIPOSOMAL DOXORUBICIN HYDROCHLORIDE", "VINCRISTINE SULFATE LIPOSOME"],
    "IO": ["TREMELIMUMAB", "BCG SOLUTION", "BCG VACCINE"],
    "Hormone": ["ANASTROZOLE", "LETROZOLE", "TAMOXIFEN CITRATE", "EXEMESTANE", "LEUPROLIDE ACETATE", "BICALUTAMIDE",
                "DEGARELIX", "MEGESTROL ACETATE", "ABIRATERONE ACETATE", "LANREOTIDE ACETATE", "OTHER HORMONE"],
    "Ancillary": ["LEUCOVORIN CALCIUM"],
    "Others": ["INVESTIGATIONAL DRUG", "OTHER ANTINEOPLASTIC", "OTHER NOS"],
}
for cls, names in CLASS_OVERRIDES.items():
    for name in names:
        drug_class[name] = cls

def get_class(raw_drug):
    base = raw_drug.split("(", 1)[0].strip().upper()   # drop the "(synonyms...)"
    if base not in drug_class:
        print(f"warning: no class for {base}, using Others")
        drug_class[base] = "Others"
    return drug_class[base]

# ---------- 2. read drug intervals (days since diagnosis) ----------
intervals = {}   # patient_id -> list of (start, end, class)
with open(DATA_PATH, newline="") as f:
    for row in csv.DictReader(f):
        if INDEX_CANCER_ONLY and row["redcap_ca_index"] != "Yes":
            continue
        for n in range(1, 6):
            drug = row[f"drugs_drug_{n}"]
            if drug == "":
                continue
            start = row[f"dx_drug_start_int_{n}"]
            if start == "":
                continue
            # end date; if blank, fall back to last administration, then to the start day
            end = row[f"dx_drug_end_int_{n}"] or row[f"dx_drug_end_or_lastadm_int_{n}"] or start
            intervals.setdefault(row["record_id"], []).append((float(start), float(end), get_class(drug)))

# ---------- 3. binary (time bin x class) features ----------
n_bins = math.ceil(HORIZON_DAYS / GRID_DAYS)
classes = sorted(set(c for drugs in intervals.values() for _, _, c in drugs))
patient_ids = list(intervals.keys())

features = np.zeros((len(patient_ids), n_bins, len(classes)), dtype=np.int8)
for p, pid in enumerate(patient_ids):
    t0 = min(s for s, _, _ in intervals[pid])            # time 0 = first treatment
    for start, end, cls in intervals[pid]:
        first_bin = int((start - t0) // GRID_DAYS)
        last_bin = int((max(end, start) - t0) // GRID_DAYS)
        if first_bin >= n_bins:                          # starts after the 5-year horizon
            continue
        features[p, first_bin:min(last_bin, n_bins - 1) + 1, classes.index(cls)] = 1

columns = [f"bin{b}_{c}" for b in range(n_bins) for c in classes]
X = pd.DataFrame(features.reshape(len(patient_ids), -1), index=patient_ids, columns=columns)
X = X.loc[:, X.var() > MIN_FEATURE_VARIANCE]
print(f"{len(patient_ids)} patients, {n_bins} bins x {len(classes)} classes, {X.shape[1]} features kept")

# ---------- 4. PCA ----------
pca = PCA(n_components=min(N_PCS, X.shape[1]), random_state=42)
scores = pca.fit_transform(X)
print(f"{pca.n_components_} PCs explain {pca.explained_variance_ratio_.sum():.0%} of the variance")

# ---------- 5. k-means for every k ----------
all_labels = pd.DataFrame(index=patient_ids)
inertia, sil = [], []
for k in K_RANGE:
    km = KMeans(n_clusters=k, n_init=50, random_state=42).fit(scores)
    all_labels[f"Cluster_k{k}"] = km.labels_ + 1
    inertia.append(km.inertia_)
    sil.append(silhouette_score(scores, km.labels_))
all_labels.to_csv("mectx_clusters.csv")

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(list(K_RANGE), inertia, "o-")
ax[0].set_title("MEC-TX k-means elbow")
ax[0].set_xlabel("k")
ax[0].set_ylabel("inertia")
ax[1].plot(list(K_RANGE), sil, "o-")
ax[1].set_title("Silhouette (higher = better)")
ax[1].set_xlabel("k")
fig.tight_layout()
fig.savefig("mectx_cluster_count.png")

# ---------- final k: plot + describe ----------
labels = all_labels[f"Cluster_k{FINAL_K}"].values

plt.figure(figsize=(7, 5))
points = plt.scatter(scores[:, 0], scores[:, 1], s=5, alpha=0.5, c=labels, cmap="tab10")
plt.legend(*points.legend_elements(), title="cluster", fontsize=7, loc="best")
plt.title(f"MEC-TX: PCA colored by k-means cluster (k={FINAL_K})")
plt.savefig("mectx_pca_clusters.png")

# heatmap per cluster: fraction of patients on each class in each 4-week bin
cluster_ids = sorted(set(labels))
fig, axes = plt.subplots(len(cluster_ids), 1, figsize=(10, 1.4 * len(cluster_ids)), sharex=True)
for ax, c in zip(axes, cluster_ids):
    frac = features[labels == c].mean(axis=0).T          # classes x bins
    ax.imshow(frac, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels(classes, fontsize=7)
    ax.set_ylabel(f"cl {c}\nn={(labels == c).sum()}", fontsize=8)
axes[-1].set_xlabel("4-week bins since first treatment")
fig.tight_layout()
fig.savefig("mectx_cluster_timelines.png")

for c in cluster_ids:
    ever = features[labels == c].max(axis=1).mean(axis=0)     # fraction ever on each class
    top = ", ".join(f"{classes[i]} {ever[i]:.0%}" for i in np.argsort(-ever)[:4])
    print(f"cluster {c} ({(labels == c).sum()} patients): {top}")

# ---------- compare with Level 1 (drug_matrix.csv k-means) if it exists ----------
try:
    level1 = pd.read_csv("drug_matrix.csv", index_col=0)
    common = [pid for pid in patient_ids if pid in level1.index]
    level1_labels = KMeans(n_clusters=FINAL_K, random_state=0, n_init=10).fit_predict(level1.loc[common])
    ari = adjusted_rand_score(level1_labels, all_labels.loc[common, f"Cluster_k{FINAL_K}"])
    print(f"agreement with Level 1 (adjusted Rand index, 0 = random, 1 = identical): {ari:.3f}")
except FileNotFoundError:
    print("drug_matrix.csv not found, skipping Level 1 comparison")
