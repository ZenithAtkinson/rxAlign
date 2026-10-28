from Patient import *
from BinnedHistory import *
class Alignments(object):
    """
    Alignments stores the alignments for each patient which were generated while getting a cluster's consensus sequence.
    While a Profile stores the counts of each event within each bin, Alignments shows how each patient's treatment path
    was aligned against all the other patient paths. 
    """
    def __init__(self):
        """
        Construct a new 'Alignments' object.
        """
        self.columns = []        # A list of dictionaries (key=Patient, value=bin) for storing alignments. Each index/column represents 1 bin_length.
        self.patients = set()    # The set of patients whose alignments are stored in self.columns

    def get_length(self):
        """
        Returns the length of the alignment in terms of bins.
        """
        return len(self.columns)
    
    def initialize_table(self, patient1:Patient, patient2:Patient, alignment1:BinnedHistory, alginment2:BinnedHistory):
        """
        Creates the inital data structure to store the first 2 alignments.

        Parameters:
        patient1 (Patient): The first patient.
        patient2 (Patient): The second patient.
        alignment1 (BinnedHistory): patient1's alignment
        alignment2 (BinnedHistory): patient2's alignment
        """
        self.patients.add(patient1)
        self.patients.add(patient2)
        for index in sorted(list(alignment1.bins.keys())):
            self.columns.append({patient1:alignment1.bins[index], patient2:alginment2.bins[index]})
    
    def update_alignment(self, profile_indicies):
        """
        Updates self.columns to in response to the new alignment between the Profile and a new BinnedHistory
        If profile_indicies contains "-" (an indel) the function will insert a new column at that index into 
        self.column

        Parameter:
        profile_indicies: a list of indicies or indels ("-") which indicate how the profile has shifted after alignment
        """
        new_length = len(profile_indicies)
        
        if new_length > self.get_length():
            for i, index in enumerate(profile_indicies):
                if index == "-":
                    self.columns.insert(i, {})
                    for patient in self.patients:
                        indel_bin = Bin()
                        indel_bin.event_names.add("-")
                        self.columns[i][patient] = indel_bin

    def add_alignment(self, patient:Patient, alignment:BinnedHistory):
        """
        Adds a patient and their aligned binned history to the alignment.

        Parameters:
        patient (Patient): patient to add to alignment
        alignment (BinnedHistory): patient's aligned binned history
        """
        self.patients.add(patient)
        
        for i, column in enumerate(self.columns):
            if len(alignment.bins[i].event_names) == 0:
                alignment.bins[i].event_names.add("-")
            column[patient] = alignment.bins[i]

    def print(self):
        """
        Prints the alignment.
        """
        for patient in self.patients:
            s = f"id: {patient.id} =>"
            for column in self.columns:
                s += column[patient].to_string()
            print(s)


        