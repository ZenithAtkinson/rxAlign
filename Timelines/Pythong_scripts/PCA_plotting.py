import os
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


df = pd.read_csv("drug_matrix.csv", index_col=0)
coords = PCA(n_components=2).fit_transform(df)

plt.figure()
plt.scatter(coords[:,0], coords[:, 1], s=5, alpha=0.4)
plt.savefig("pca.png")

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
fig.savefig("cluster_count.png")

# --- final clusters: k-means ---
# (hierarchical + jaccard put ~1390 patients in one cluster and split off single
#  patients with rare drugs, so its silhouette scores above are misleading)
K = 6
labels = KMeans(n_clusters=K, random_state=0, n_init=10).fit_predict(df) + 1

plt.figure(figsize=(7, 5))
points = plt.scatter(coords[:, 0], coords[:, 1], s=5, alpha=0.5, c=labels, cmap="tab10")
plt.legend(*points.legend_elements(), title="cluster", fontsize=7, loc="best")
plt.title(f"PCA colored by k-means cluster (k={K})")
plt.savefig("pca_clusters.png")

# what's in each cluster: size + drugs given to the most patients in it
short_names = [c.split("(")[0] for c in df.columns]   # drop the "(synonyms...)" for printing
for c in sorted(set(labels)):
    in_cluster = df[labels == c]
    share = pd.Series(in_cluster.mean().values, index=short_names).sort_values(ascending=False)
    top = ", ".join(f"{name} {frac:.0%}" for name, frac in share.head(4).items())
    print(f"cluster {c} ({len(in_cluster)} patients): {top}")

#umap plot? with umpa_learn

#check cluster mean something after with the (os_d_status) between clusters. Sad but compares final results, health, and survival

#THEN treatment order?
