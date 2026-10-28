from Cluster import *
from Patient import *
import numpy as np


class DistanceMatrix(object):
    """
    Data structure for storing patient to patient distances.
    """
    def __init__(self, patients):
        """
        Construct a new 'DistanceMatrix' object. 

        Parameter:
        patients: set of 'Patient' objects from which 'DistanceMatrix' object will be constructed
        """
        self.patients = patients
        self.patient_to_index_map = {}
        self.index_to_patient_map = {}
        self.matrix = []

        # initialize the matrix structure
        for i in range(len(patients)):
            self.matrix.append([0]*len(patients))
        
        self.index_patients()

    def get_distance(self, patient1:Patient, patient2:Patient):
        """
        Returns the distance between two patients.
        """
        i = self.patient_to_index_map[patient1]
        j = self.patient_to_index_map[patient2]
        if type(i) == str:
            i = int(i)
            j = int(j)
        return self.matrix[i][j]   

    def index_patients(self):
        """
        Loops through self.patients  and assigns indices to each patient. The
        indices are stored in self.patient_to_index_map and self.index_to_patient_map
        to enable easy look-up.
        """
        for i, patient in enumerate(self.patients):
            self.patient_to_index_map[patient] = i
            self.index_to_patient_map[i] = patient

    def add_row(self, patient:Patient, distance_dict:dict):
        """
        Adds a patient and their distance scores to matrix.
        
        Parameters:
        patient (Patient): patient to add to 'DistanceMatrix'
        distance_dict (dict): dict of distances between patient and other patients
        in matrix (key = other_patient (Patient), value = distance)
        """
        self.patients.add(patient) # may need to delete
        i = len(self.patients)-1
        self.patient_to_index_map[patient] = i
        self.index_to_patient_map[i] = patient
        row = [0]*(len(self.patients))
        # self.matrix = np.vstack([self.matrix, row])
        self.matrix.append(row)
        for patient2 in self.patients:
            j = self.patient_to_index_map[patient2]
            
            self.matrix[i][j] = distance_dict[patient2]
            if i != j:
                self.matrix[j].append(distance_dict[patient2])
        
    def record_distance(self, patient1:Patient, patient2:Patient, distance):
        """
        Records the distance between patient1 and patient2 in matrix.
        """
        i = self.patient_to_index_map[patient1]
        j = self.patient_to_index_map[patient2]
        self.matrix[i][j] = distance
        self.matrix[j][i] = distance

    def write(self, file_path):
        """
        Writes distance matrix to a tsv.

        Parameter:
        file_path: file_path (tsv) of where to save distance matrix
        """
        cm = self.matrix
        with open(file_path, 'w') as output:
            header = ""
            for i in range(len(cm)):
                if i not in self.index_to_patient_map.keys():
                    i = str(i)
                header += f"{self.index_to_patient_map[i].id}\t"
            header += '\n'
            output.write(header)
            for i, row in enumerate(cm):
                s = ""
                if i not in self.index_to_patient_map.keys():
                    i = str(i)
                s += f"{self.index_to_patient_map[i].id}\t"
                for j, distance in enumerate(row):
                    s += f"{round(distance, 3)}\t"
                output.write(s)
                output.write('\n')
    
    def filter_by_cluster(self, cluster, new_patient):
        """
        Returns a new 'DistanceMatrix' which only contains distances between patients
        from the specified cluster
        """
        filtered_matrix = DistanceMatrix(cluster.patients)

        for patient1 in cluster.patients:
            for patient2 in cluster.patients:
                if patient1 != new_patient and patient2 != new_patient:
                # TODO: could make faster
                    i = self.patient_to_index_map[patient1]
                    j = self.patient_to_index_map[patient2]
                    if type(i) == str:
                        i = int(i)
                        j = int(j)
                    distance = self.matrix[i][j]
                    filtered_matrix.record_distance(patient1, patient2, distance)
            
        return filtered_matrix


