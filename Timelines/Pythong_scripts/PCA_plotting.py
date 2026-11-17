import os
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#using python for sklearn and scipy

# folder this script is in (Pythong_scripts), so files are found/saved there no matter where it's run from
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

MATRIX_FILE = "newer_timeline_matrix.csv"

MATRIX_FILE = "newer_timeline_matrix.csv"   # Level 2; use "drug_matrix.csv" for Level 1
NAME = MATRIX_FILE.replace(".csv", "")      # goes in the plot filenames so runs don't overwrite each other

df = pd.read_csv(os.path.join(SCRIPT_DIR, MATRIX_FILE), index_col=0)
coords = PCA(n_components=2).fit_transform(df)

plt.figure()
plt.scatter(coords[:,0], coords[:, 1], s=5, alpha=0.4)
plt.savefig(os.path.join(SCRIPT_DIR, f"{NAME}_pca.png"))

ks = range(2, 16)

# --- k-means: elbow plot + silhouette ---
inertia = []
kmeans_sil = []
for k in ks:
    km = KMeans(n_clusters=k, random_state=0, n_init=10).fit(df)
    inertia.append(km.inertia_)
    kmeans_sil.append(silhouette_score(df, km.labels_))

# --- hierarchical (jaccard distance, average linkage): silhouette ---
dist = pdist(df.values.astype(bool), metric="jaccard")   # distance between every pair of patients
D = squareform(dist)
Z = linkage(dist, method="average")                      # build the cluster tree
hier_sil = []
for k in ks:
    labels = fcluster(Z, t=k, criterion="maxclust")
    # the top merges are tied at distance 1.0 (no shared drugs), so some k can't be cut exactly
    if len(set(labels)) < 2:
        hier_sil.append(np.nan)
        continue
    hier_sil.append(silhouette_score(D, labels, metric="precomputed"))

print("k\tkmeans_inertia\tkmeans_silhouette\thierarchical_silhouette")
for k, i, s1, s2 in zip(ks, inertia, kmeans_sil, hier_sil):
    print(f"{k}\t{i:.1f}\t{s1:.3f}\t{s2:.3f}")

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(ks, inertia, "o-")
ax[0].set_title("K-means elbow")
ax[0].set_xlabel("k")
ax[0].set_ylabel("inertia")
ax[1].plot(ks, kmeans_sil, "o-", label="k-means")
ax[1].plot(ks, hier_sil, "o-", label="hierarchical (jaccard)")
ax[1].set_title("Silhouette (higher = better)")
ax[1].set_xlabel("k")
ax[1].legend()
fig.tight_layout()
fig.savefig(os.path.join(SCRIPT_DIR, f"{NAME}_cluster_count.png"))

# --- final clusters: k-means ---
# (hierarchical + jaccard put ~1390 patients in one cluster and split off single
#  patients with rare drugs, so its silhouette scores above are misleading)
# i dont see a clear elbow, maybe 5? k-means silhouette is highest at 15 (Level 1); 11 chosen as readable.
K = 11
labels = KMeans(n_clusters=K, random_state=0, n_init=10).fit_predict(df) + 1

plt.figure(figsize=(7, 5))
points = plt.scatter(coords[:, 0], coords[:, 1], s=5, alpha=0.5, c=labels, cmap="tab10")
plt.legend(*points.legend_elements(), title="cluster", fontsize=7, loc="best")
plt.title(f"PCA colored by k-means cluster (k={K})")
plt.savefig(os.path.join(SCRIPT_DIR, f"{NAME}_pca_clusters.png"))

# what's in each cluster: size + drugs given to the most patients in it
short_names = [c.split("(")[0] for c in df.columns]   # drop the "(synonyms...)" for printing
for c in sorted(set(labels)):
    in_cluster = df[labels == c]
    share = pd.Series(in_cluster.mean().values, index=short_names).sort_values(ascending=False)
    top = ", ".join(f"{name} {frac:.0%}" for name, frac in share.head(4).items())
    print(f"cluster {c} ({len(in_cluster)} patients): {top}")

# --- agreement with MEC-TX (adjusted Rand index: 0 = chance, 1 = identical groupings) ---
# Shows how much these drug-based clusters match MEC-TX's time-based ones, so the number in the
# report can be re-run. Needs R_scripts/mectx_clusters.csv (run mectx_run.R first); compares
# patients in both files (MEC-TX drops some) against MEC-TX at the same k.
from sklearn.metrics import adjusted_rand_score
mectx = pd.read_csv(os.path.join(SCRIPT_DIR, "..", "R_scripts", "mectx_clusters.csv"), index_col="sample")[f"Cluster_k{K}"]
level1 = pd.Series(labels, index=df.index)
common = level1.index.intersection(mectx.index)
print(f"ARI vs MEC-TX (k={K}, {len(common)} patients): {adjusted_rand_score(level1[common], mectx[common]):.3f}")

#umap plot? with umpa_learn

#check cluster mean something after with the (os_d_status) between clusters. Sad but compares final results, health, and survival

#THEN treatment order?
