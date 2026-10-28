from DistanceMatrix import *
from Profile import * 

class Cluster(object):
    """
    Stores similar patients and their alignment.
    """
    def __init__(self, id:int):
        """
        Construct a 'Cluster' object

        Parameter:
        id (int): cluster id
        """
        self.id = id
        self.patients = set()       # set of patients that belong to this cluster
        self.profile = Profile()    # profile for storing the patient's alignment. 
        self.dist_data = []         # transformed distance matrix for patients of the cluster for easily getting the next most similar patient 
    
    def add_patient(self, patient):
        """
        Adds a patient to self.patients.
        """
        self.patients.add(patient)

    def fill_in_distances(self, original_distance_matrix):
        """
        Populates self.dist_data with the distances from original_distance_matrix for the patients within this cluster 
        
        Parameter:
        original_distance_matrix (DistanceMatrix): DistanceMatrix with distances for all patients
        """
        self.dist_data = []
        for patient1 in self.patients:

            for patient2 in self.patients:
                # this is twice as big as it needs to be but the code is easier if we keep it like this
                if patient1 != patient2:

                    distance = original_distance_matrix.get_distance(patient1, patient2)                    
                    label = self.make_unique_patient_to_patient_label(patient1_id=patient1.id, patient2_id=patient2.id)
                    self.dist_data.append({'patient1':patient1, 'patient2':patient2, 'distance':distance, 'label':label})
                
    def make_unique_patient_to_patient_label(self, patient1_id, patient2_id):
        """
        Returns a unique label for a patient to patient comparison
        """
        ids = sorted([patient1_id, patient2_id])
        label = f"{ids[0]}to{ids[1]}"
        return label
    
    def get_starting_patients(self, num_starts:int):
        """
        Returns an ordered list of the num_starts-most similar pairs of patients

        Parameter:
        num_starts(int): number of pairs of patients to return
        """
        self.dist_data = sorted(self.dist_data, key=lambda d: d['distance'], reverse=False)
        pairs = []
        labels = set()      # set of patient comparisons already in pairs. Used to make sure not duplicate pairs are returned.
        length = num_starts

        # adjust number of pairs to return if nessesary
        if len(self.dist_data) < num_starts:
            length = len(self.dist_data)
            print(f"new length {length}")

        i = 0
        while len(pairs) < length:
            if self.dist_data[i]['label'] not in labels:
                pair = [self.dist_data[i]['patient1'], self.dist_data[i]['patient2']]
                pairs.append(pair)
                labels.add(self.dist_data[i]['label'])
            i += 1

        return pairs