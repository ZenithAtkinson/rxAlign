import numpy as np

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform
from DistanceMatrix import *
from Cluster import *

import matplotlib.pyplot as plt

class HierarchicalClusterer:
    """
    Clusters patients and plots dendrogram
    """
    def get_cluster_labels(self, labels):
        """
        Returns a list of labels that contains each cluster name only once.

        Parameter:
        labels (list): list of labels returned from fcluster()
        """
        seen_clusters = set()
        cluster_labels = []
        for i, label in enumerate(labels):
            if label not in seen_clusters:
                seen_clusters.add(label)
                cluster_labels.append(label)
            else:
                cluster_labels.append(None)
        
        return cluster_labels
    
    def get_patient_labels(self, distance_matrix:DistanceMatrix, query_patient=None):
        """
        Returns a list of patient ids for the dendrogram.
        """
        labels = []
        for i in sorted(distance_matrix.index_to_patient_map.keys()):
            patient_id = distance_matrix.index_to_patient_map[i].id
            if patient_id == query_patient:
                patient_id = "Query: " + patient_id
            labels.append(patient_id)
        return labels

    def cluster(self, distance_matrix:DistanceMatrix, linkage_method, distance_threshold, return_Z=False):
        """
        Clusters patients using hierarchical clustering.

        Parameters:
        distance_matrix (DistanceMatrix): distance matrix containing patient to patient distance distances
        linkage_method: "ward", "complete", "single", "average", or "centroid"
        distance_threshold: threshold above which clusters will not be formed
        return_Z: return model
        """
        dists = squareform(distance_matrix.matrix)
        Z = linkage(dists, linkage_method)
        cluster_labels = fcluster(Z, t=distance_threshold, criterion='distance')
        
        if return_Z:
            return cluster_labels, Z
        else:
            return cluster_labels
        
    def cluster_and_plot(self, distance_matrix:DistanceMatrix, linkage_method, distance_threshold, title, file_path, show_id=False, query_patient=False):
        """
        Clusters patients and plots dendrogram.

        Parameters:
        distance_matrix (DistanceMatrix): distance matrix containing patient to patient distance distances
        linkage_method: "ward", "complete", "single", "average", or "centroid"
        distance_threshold: threshold above which clusters will not be formed
        title: dendrogram title
        file_path: where to save dendrogram
        show_id: if True, dendrogram will show patient ids. If False, will show cluster labels
        """
        cluster_labels, Z= self.cluster(distance_matrix=distance_matrix, linkage_method=linkage_method, distance_threshold=distance_threshold, return_Z=True)
        plt.figure()  # just a way to make sure that we start with an empty figure

        if show_id:
            dendrogram(Z=Z, labels=self.get_patient_labels(distance_matrix=distance_matrix, query_patient=query_patient),above_threshold_color="black", color_threshold=distance_threshold)
            plt.xticks(rotation=90)
        else:
            dendrogram(Z=Z, labels=self.get_cluster_labels(labels=cluster_labels),above_threshold_color="black", color_threshold=distance_threshold)
        
        
        plt.title(title)
        
        if file_path:
            plt.savefig(file_path, bbox_inches = "tight")
        
        
        return cluster_labels

