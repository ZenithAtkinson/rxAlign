import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
import time
import copy
import os.path




from DistanceMatrix import *
from LongestCommonSubsequenceSolver import *
from Cluster import *
from FileReader import *
from HierarchicalClusterer import *
from Encoder import *

class Controller:
    def __init__(self):
        self.clusters = None
        self.lcs = None

    def run_all(self, 
                patient_data_path,     
                med_to_med_class_path,
                bin_length=None,   
                use_treatment_lines=False, 
                regimen_column=None,
                combo_cutoff=0.01,   
                MATCH_REWARD=-2,       
                MISMATCH_PENALTY=1,
                INDEL_PENALTY=0.1, 
                merging_clusters_dict={},
                consensus_seq_starts=2,
                consensus_seq_sub_starts=5,
                linkage_type='ward',
                distance_threshold=2,   # distance threshold above which a cluster will not be formed
                cs_MATCH_REWARD=-10,
                cs_MISMATCH_PENALTY=0,
                cs_INDEL_PENALTY=4,
                seq_logo_cut_off=0.05,  # frequency threshold below which a med will not be shown
                med_letters=4,
                save_data=False,
                folder=None):
        
        """
        Reads in files, bins patient histories, aligns patient histories, clusters patients and plots dendrogram, generates and plots
        consensus sequences for each cluster.

        Parameters:
        patient_data_path:  path to file with patient history data
        med_to_med_class_path: path to file with medication to med class
        bin_length: how long bins should be (unit is the same as unit of start and end dates in input file)
        use_treatment_lines: if True, uses treatment lines instead of time binning. Specify regimen column (or treatment line column) if desired.
        regimen_column: Name of column in input data that specifies treatment line. Only used if 'use_treatment_lines' is True.
        combo_cutoff: Cutoff (in units of start and end dates from input data) for determining combination therapy. Used if regimen_column is None.
        MATCH_REWARD (int or float): Reward for matching events. Should be negative.
        MISMATCH_PENALTY (int or float): Penalty for an event not occuring in both bins. Should be positive.
        INDEL_PENALTY (int or float): -INDEL_PENALTY is used as the distance of two empty bins. Should be positive.    
        merging_clusters_dict (dict): dictionary of cluster_ids to cluster_ids. Enables user to merge clusters generated from clustering
        consensus_seq_starts (int): number of starts to use when generating consensus sequences
        consensus_seq_sub_starts (int): number of substarts to use when generating consensus sequences
        linkage_type='ward': linkage type to use in clustering model (must be "ward", "complete", "single", "average", or "centroid")
        distance_threshold (float or int):  distance threshold above which a cluster will not be formed
        seq_logo_cut_off (float): frequency threshold below which a med will not be shown
        save_data (bool): if True, saves data as json
        folder: folder of where to save the data. Create folder before running.
        """
        fr = FileReader()
        self.color_dict = fr.get_color_dict(med_to_med_class_path, length=med_letters)

        self.patients = fr.get_binned_patient_histories(
            file_path=patient_data_path, 
            bin_length=bin_length, 
            use_treatment_lines=use_treatment_lines,
            regimen_column=regimen_column,
            combo_cutoff=combo_cutoff)
        
        
        self.distance_matrix = self.align_patients(
            self.patients, 
            MATCH_REWARD=MATCH_REWARD,
            MISMATCH_PENALTY=MISMATCH_PENALTY,
            INDEL_PENALTY=INDEL_PENALTY
            )
        
        self.distance_matrix.write((folder+ r'\distance_matrix.tsv')
                                   )
        if folder:
           dendrogram_file_path= folder + r'\dendrogram.png'
        else:
           dendrogram_file_path = 'dendrogram.png'

        self.clusters = self.cluster_patients(
            distance_matrix=self.distance_matrix,
            linkage_type=linkage_type,
            distance_threshold=distance_threshold,
            file_path=dendrogram_file_path,
            merging_clusters_dict=merging_clusters_dict
            )
        
        if folder:
            cluster_file_path = folder + r"\data_with_clusters.tsv"
        else:
            cluster_file_path = r"\data_with_clusters.tsv"

        self.id_to_patient_dict = {}
        for patient in self.patients:
            self.id_to_patient_dict[patient.id] = patient
        og_data = pd.read_csv(patient_data_path, sep='\t')
        og_data["Cluster"] = None
        for index, row in og_data.iterrows():
            patient = self.id_to_patient_dict[row["ID"]]              
            og_data.at[index, 'Cluster'] = patient.cluster
        og_data.to_csv(cluster_file_path, sep='\t')

        lcs = LongestCommonSubsequenceSolver(
            MATCH_REWARD=cs_MATCH_REWARD,
            MISMATCH_PENALTY=cs_MISMATCH_PENALTY,
            INDEL_PENALTY=cs_INDEL_PENALTY
        )
        for id in sorted(self.clusters.keys()):
            print(id)
            if id in self.clusters.keys():
                cluster = self.clusters[id]
                profile = self.get_consensus_sequences(
                    cluster, 
                    num_starts=consensus_seq_starts, 
                    num_sub_tries=consensus_seq_sub_starts,
                    lcs=lcs
                    )
                
                print(profile.get_simplified_string(threshold=seq_logo_cut_off))
                
                seq_logo_title = f"cluster_{cluster.id}"
                if folder:
                    seq_logo_file_path = folder + r"\seq_logo" + seq_logo_title +".png"
                else:
                    seq_logo_file_path = "seq_logo" + seq_logo_title +".png"

                profile.get_logo(file_path=seq_logo_file_path, 
                                color_dict=self.color_dict, 
                                threshold=seq_logo_cut_off,
                                title=seq_logo_title,
                                med_letters=med_letters)
        
        if save_data:
            self.save_data(folder_path=r"/instance") # where to save data
            
         
        
    def align_patients(self, patients, MATCH_REWARD=-2, MISMATCH_PENALTY=1,INDEL_PENALTY=0.1):
        # somehow the indel penalty needs to decrease with decreasing bin length
        # take avg history length and divide it by bin length...
        self.patients = patients
        self.distance_matrix = DistanceMatrix(patients)
        self.lcs = LongestCommonSubsequenceSolver(MATCH_REWARD=MATCH_REWARD, MISMATCH_PENALTY=MISMATCH_PENALTY,INDEL_PENALTY=INDEL_PENALTY)
        total = len(patients)**2
        complete = 0
        for i, patient1 in enumerate(patients):
            for j, patient2 in enumerate(patients[i+1:]):
                distance = self.lcs.getDistanceandAlignment(patient1.binned_history, patient2.binned_history)
                self.distance_matrix.record_distance(patient1, patient2, distance) # fix, change to use id
                # complete += 1
                # if complete/total
                # print(complete/total * 100)
        return self.distance_matrix

    def cluster_patients(self, distance_matrix:DistanceMatrix, linkage_type, file_path=False, title="All Clusters", distance_threshold=0.5, merging_clusters_dict={}):
        # prepare matrix
        # TODO: need to figure out what the distance threshold number needs to be
        
        # ap = AgglomerativeClustering(n_clusters=None, metric='precomputed', linkage=linkage_type, distance_threshold=distance_threshold)
        # self.model = ap # TODO: make this cleaner by refactoring
        # # Fit the model to the distance matrix
        # ap.fit(distance_matrix.matrix)
        
        # cluster_labels = ap.labels_
        # num_clusters_ap = ap.n_clusters_
    

        # # num_clusters_ap = len(ap.cluster_centers_indices_)
        # print("Number of clusters (Affinity Propagation):", num_clusters_ap)
        # print("\nAffinity Propagation Cluster Assignments:")
        
        clusterer = HierarchicalClusterer()
        cluster_labels = clusterer.cluster_and_plot(
            distance_matrix=distance_matrix, 
            linkage_method=linkage_type, 
            distance_threshold=distance_threshold, 
            title=title, 
            file_path=file_path)

        clusters = {}
        
        for i, cluster in enumerate(cluster_labels):
            if i not in distance_matrix.index_to_patient_map.keys():
                i = str(i)
            patient = distance_matrix.index_to_patient_map[i]

            if cluster in merging_clusters_dict.keys():
                cluster = merging_clusters_dict[cluster]

            if cluster not in clusters.keys():
                clusters[cluster] = Cluster(cluster)

            patient.cluster = cluster
            clusters[cluster].add_patient(patient)
        self.clusters = clusters
        print(len(clusters))
        return clusters
        
    def query(self, 
              query_file, 
              dendrogram_path, 
              med_to_med_class_path, 
              use_treatment_lines=True, 
              bin_length=False, 
              regimen_column=False, 
              combo_cutoff=False, 
              linkage_type="ward"):

        fr = FileReader()
        self.color_dict = fr.get_color_dict(med_to_med_class_path, length=4)

        if self.lcs == None:
            self.lcs = LongestCommonSubsequenceSolver()

        patient  = fr.get_binned_patient_histories(
            file_path=query_file, 
            bin_length=bin_length, 
            use_treatment_lines=use_treatment_lines,
            regimen_column=regimen_column,
            combo_cutoff=combo_cutoff)[0]
        
        min_dist = float("inf")
        matching_cluster = None
        for id, cluster in self.clusters.items():
            alignment, profile_indices, distance = self.lcs.alignWithProfile(patient.binned_history, cluster.profile)
            # print(f'cluster:{cluster.id}, distance:{distance}')
            # cluster.profile.alignments.print()
            # cluster.profile.update_alignment(profile_indices)   # *new_profile.num_histories # does it make sense to do this?
            # cluster.profile.add_history(alignment, patient)
            # cluster.profile.alignments.print()
            # alignment.to_string()
            if distance < min_dist:
                min_dist = distance
                matching_cluster = cluster
        
        distance_dict = {}
        matching_cluster.patients.add(patient)
        print(f'cluster:{matching_cluster.id}, distance:{min_dist}')
        for patient1 in matching_cluster.patients:
            distance = self.lcs.getDistanceandAlignment(patient.binned_history, patient1.binned_history)
            distance_dict[patient1] = distance
        
        new_matrix = self.distance_matrix.filter_by_cluster(matching_cluster, patient)
        for p in distance_dict.keys():
            new_matrix.record_distance(patient, p, distance_dict[p])
        
        clusterer = HierarchicalClusterer()
        title = f"{patient.id} in cluster {matching_cluster.id}"
        clusterer.cluster_and_plot(distance_matrix=new_matrix, linkage_method=linkage_type, distance_threshold=0, title=title, file_path=dendrogram_path, show_id=True, query_patient=patient.id)
        
        
        print(matching_cluster.id)
        matching_cluster.patients.remove(patient)

    def get_consensus_sequences(self, cluster:Cluster, num_starts=5, num_sub_tries=5, lcs=None):
        if lcs != None:
            self.lcs = lcs
        # MATCH_REWARD=-2, MISMATCH_PENALTY=1,INDEL_PENALTY=0.1
        if len(cluster.patients) == 0:
            print("outlier")
            # need a funtion that changes a history into a profile/alignment
        updates = 0

        if len(cluster.patients) == 1:
            return False
        cluster.fill_in_distances(self.distance_matrix)
        total_distance = 0 
        # get a top num_starts list of patients to try to create the matrix with

        pairs = cluster.get_starting_patients(num_starts)
        best_distance = float("inf")
        best_profile = None
        for pair in pairs:
            new_profile = Profile(cluster.dist_data)
            patient1 = pair[0]
            patient2 = pair[1]
            alignment1, alignment2, distance = self.lcs.getDistanceandAlignment(patient1.binned_history,patient2.binned_history,returnAlignment=True)
            new_profile.add_history(alignment1, patient1)
            new_profile.add_history(alignment2, patient2)
            new_profile.alignments.initialize_table(patient1, patient2, alignment1, alignment2)
            new_profile.distance += distance
            # print("Starting patients:")
            # print(patient1.id)
            # print(patient1.binned_history.to_string())
            # print(patient2.id)
            # print(patient2.binned_history.to_string())
            # print("alignment")
            # print(new_profile.get_simplified_string(threshold=0))
            # new_profile.alignments.print()
            # print("profile")
            # new_profile.print()
            
            # print(f"distance: {distance}")


            while new_profile.num_histories < len(cluster.patients):
                next_patients = new_profile.get_next_patient(num_sub_tries)
                best_sub_distance = float("inf")
                best_alignment = {} # should be [alignment, profile, patient]
                for next_patient in next_patients:
                    # print(f"next similar patient: {next_patient.id}")
                    # print(next_patient.binned_history.to_string())

                    # TODO: make sure this isn't messing anything up within the function 
                    alignment, profile_indices, sub_distance = self.lcs.alignWithProfile(next_patient.binned_history, new_profile)
                    # print("aligned history:")
                    # print(alignment.to_string())
                    # print(f"distance: {sub_distance}")
                    if sub_distance < best_sub_distance:
                        # print(f"found better sub alignment: old distance: {best_sub_distance}, new_distance: {sub_distance}")
                        best_alignment = {'alignment':alignment, 'profile_indices':profile_indices, 'next_patient':next_patient}
                        best_sub_distance = sub_distance
                    
                new_profile.distance += best_sub_distance
                if new_profile.distance > best_distance: # kind of like TSP => this is like TSP, isn't? so does that mean that this approach is flawed because I'm not making it random enough if choosing which state to go to?
                    # print(f"distance {new_profile.distance} too high, aborting iteration")
                    break # TODO: check logic here
                else:
                    # print(f"finished iterations, adding patient {best_alignment['next_patient'].id} and distance is now {new_profile.distance}")
                    new_profile.update_alignment(best_alignment['profile_indices'])   # *new_profile.num_histories # does it make sense to do this?
                    new_profile.add_history(best_alignment['alignment'], best_alignment['next_patient'])
                    # print("alignment")
                    # new_profile.alignments.print()
                    # print(new_profile.get_simplified_string(threshold=0))
                    # print("profile")
                    # new_profile.print()
            
            if new_profile.distance < best_distance:
                # print(f"found better alignment: old distance: {best_distance}, new_distance: {new_profile.distance}")
                updates += 1
                best_profile = new_profile
                best_distance = new_profile.distance
            
        cluster.profile = best_profile
        # add the best profile to the cluster
        # print("best profile")
        # cluster.profile.alignments.print()
        # print(f"distance: {best_distance}")
        # print(f"updates: {updates}")
        return best_profile



        # need to get the next patient to align
        # store the sorted similarity distance in a sorted something then iterate through the list, starting with the 
        # most similar distance, check if patients are in the cluster and if one of the patients is in the alignment

        # need to create an alignment algorithm that works with a profile
    def save_data(self, folder_path, bin_length=False):
        # things to save: distance matrix, patients, clusters (how to we make sure we're not saving the patient data multiple times?)
        encoder = Encoder()
        if bin_length:
            encoder.save(self.distance_matrix, (folder_path + r'\distance_matrix_'+str(bin_length) + '.json'))
        else:
            encoder.save(self.distance_matrix, (folder_path + r'\distance_matrix.json'))
        if self.clusters != None:
            for id, cluster in self.clusters.items():
                cluster.profile.patients = None
                cluster.profile.alignments.patients = None
                cluster.dist_data = None
                cluster.profile.dist_data = None
                # switch out patient objects for patient ids in alignments columns
                columns = []
                for column in cluster.profile.alignments.columns:
                    new_column = {}
                    for patient, bin in column.items():
                        new_column[patient.id] = bin
                    columns.append(new_column)
                cluster.profile.alignments.columns = columns
            
            encoder.save(self.clusters, (folder_path + r'\clusters.json'))
        
    def load_data(self, folder_path, bin_length=False):
        encoder = Encoder()
        # check if file exists in folder first
        if bin_length:
            self.distance_matrix = encoder.load((folder_path + r'\distance_matrix_'+str(bin_length) + '.json'))
        else:
            self.distance_matrix = encoder.load((folder_path + r'\distance_matrix.json'))
        self.patients = copy.copy(self.distance_matrix.patients)
        # fix indicies to be int intead of str TODO: find a better solution
        for p in self.patients:
            new_bins = {}
            for index, bin in p.binned_history.bins.items():
                new_bins[int(index)] = bin
            p.binned_history.bins = new_bins
        id_to_patient_dict = {}
        for patient in self.patients:
            id_to_patient_dict[patient.id] = patient
        self.id_to_patient_dict = id_to_patient_dict
        new_patient_to_index_map = {}
        for index, patient in self.distance_matrix.index_to_patient_map.items():
            patient_object = id_to_patient_dict[patient.id]
            self.distance_matrix.index_to_patient_map[index] = patient_object
            new_patient_to_index_map[patient_object] = index
        self.distance_matrix.patient_to_index_map = new_patient_to_index_map

        if os.path.exists((folder_path + r'\clusters.json')):
            self.clusters = encoder.load((folder_path + r'\clusters.json'))
       
            for id, cluster in self.clusters.items():
                new_set = set()
                for patient in cluster.patients:
                    new_set.add(id_to_patient_dict[patient.id])
                cluster.patients = new_set
                
                cluster.profile.patients = copy.copy(cluster.patients)
                cluster.profile.alignments.patients = copy.copy(cluster.patients)
                # pick back up here later
                columns = []
                for column in cluster.profile.alignments.columns:
                    new_column = {}
                    for id, bin in column.items():
                        new_column[id_to_patient_dict[id]] = bin
                    columns.append(new_column)
                cluster.profile.alignments.columns = columns
            
        

        

if __name__ == '__main__':
    

    start_time = time.time()
    # Your main code here
    
    controller = Controller()
    fr = FileReader()
    patients = fr.get_binned_patient_histories('data/MEDICATION_DATA_FINAL.txt', 0.25) 
    distance_matrix = controller.align_patients(patients)
    controller.cluster_patients(distance_matrix)

    with open('data/PatientCluster.tsv', 'w') as output:
        output.write("PATIENTID\tCluster\n")
        for patient in patients:
            output.write(f"{patient.id}\t{patient.cluster}\n")
    
    elapsed_time = time.time() - start_time

    print(f"Execution time: {elapsed_time:.2f} seconds")

# now need to get this to work with the sample data and try it out!
