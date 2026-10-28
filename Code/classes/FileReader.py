from Patient import *
from Event import *
import pandas as pd
import math
import copy

class FileReader:
    """
    Reads file into patient objects with binned histories
    """

    def __init__(self):
        self.seen_meds = set()
    
    def get_med_abbrv(self, med_name, length):
        if len(med_name) >= 3:
            abbrv = med_name[:length].upper()
            
        else:
            abbrv = copy.copy(med_name)
            while len(abbrv) < length:
                abbrv += "-"
        
        return abbrv


    def get_color_dict(self, file_path, length):
        self.abbrv_length = length
        med_class_info = pd.read_csv(file_path, sep='\t')
        class_color_dict = {
            'Targeted': 'blue',
            'Chemo': 'red',
            'ICB': 'green',
            'Other': 'gray'
        }
        self.color_dict = {}
        for index, row in med_class_info.iterrows():
            medication = row['Medication'].upper()
            self.seen_meds.add(medication)
            abbrv = self.get_med_abbrv(medication, length)
            self.color_dict[abbrv] = class_color_dict[row['Class']]
        return self.color_dict
    
    def get_binned_patient_histories(self, file_path, bin_length=None, num_bins=5, use_treatment_lines=False, combo_cutoff=0.01, regimen_column=None):
        """
        Creates a list of patient objects from file_path

        Parameters:
        file_path: path to patient data file
        bin_length (int or float): length of bins patient histories will be binned into
        num_bins (int): if bin_length is not provided, will bin histories into specified number of bins
        use_treatment_lines (bool): run algorithm using treatment lines rather than binning
        combo_cutoff (int or float): cutoff below which two medications are considered part of the same treatment line
        """
        data = pd.read_csv(file_path, sep='\t')
        data = data.sort_values(by=['ID', 'Start'])
        patients = []

        missing_meds_from_med_class = set()

        if not use_treatment_lines and bin_length is None:
            bin_length = self.calculate_bin_length(data, num_bins)

        for i, id in enumerate(data['ID'].unique()):

            patient_data = data.loc[data['ID'] == id]

            if len(patient_data) == 0:
                continue

            age_at_treatment_start = float(patient_data.iloc[0]['Start'])
            if use_treatment_lines:
                bin_length = 1
            new_patient = Patient(id, bin_length)
            patients.append(new_patient)

            last_bin = -1 # or last treament line
            last_start = - (combo_cutoff)*2
            for index, row in patient_data.iterrows():
                start = float(row['Start'])-age_at_treatment_start
                end = row['End']-age_at_treatment_start
                if type(end) is not int and type(end) is not float:
                    if index < (len(patient_data)-1):
                        end = patient_data.iloc[index+1]['Start']-1
                        if end < start:
                            end = start+0.00001
                    else:
                        end = start+0.00001
                
                med_data = row['Medication']
                if "," in row["Medication"]:
                    events = med_data.split(", ")
                
                else:
                    events = [med_data]
                
                for event in events:
                    event = event.upper()
                    if event not in self.seen_meds and self.get_med_abbrv(event, self.abbrv_length) not in self.color_dict.keys():
                        missing_meds_from_med_class.add(event)
                    if use_treatment_lines:
                        if (start - last_start) >= combo_cutoff:
                            last_bin += 1
                            last_start = start

                        new_event = Event(
                            name=event,
                            start=last_bin,
                            end=last_bin
                            )
                        
                        if regimen_column:
                            bin = row[regimen_column]-1
                            
                            if str(bin) != "nan":                      
                                new_event = Event(
                                name=event,
                                start=bin,
                                end=bin
                                )
                                           
                    else:                             
                        # use time binning
                        new_event = Event(
                            name=event,
                            start=start,
                            end=start
                            )
                    new_patient.add_event(new_event)

        for med in missing_meds_from_med_class:

            print(f"error: {med} not in med class file, setting color to black")
            self.color_dict[self.get_med_abbrv(med_name=med, length=self.abbrv_length)] = "black"

        return patients

    def calculate_bin_length(self, data, num_bins=5):

        lengths = []
        for i, id in enumerate(data['ID'].unique()):
            patient_data = data.loc[data['ID'] == id]
            age_at_treatment_start = float(patient_data.iloc[[0]]['Start'])

            age_at_treatment_end = float(patient_data.iloc[[-1]]['End'])
            if math.isnan(age_at_treatment_end):
                age_at_treatment_end = age_at_treatment_start+0.0001
            length = age_at_treatment_end - age_at_treatment_start
            lengths.append(length)
        avg_length = sum(lengths)/len(lengths)
        bin_length = avg_length//num_bins
         
        return bin_length

    