from Patient import *
from Bin import *
from BinnedHistory import *
from Alignments import *
from Med_Logo_Maker.Logo import *
import copy


class Profile(object):
    """
    Data structure for storing the profile and alignments of a cluster's patients.
    """
    def __init__(self, dist_data=None):
        self.length = 0
        self.num_histories = 0  # number of histories that have been added to the profile
        self.bin_columns = []   # list of 'Bins' in profile
        self.patients = set()   # set of patients whose histories have been added to the profile
        self.distance = 0       # the sum of all the distances returned from aligning the profile to the next history added in 
        self.alignments = Alignments() # stores each history and how it is aligned with the other histories
        if dist_data:
            self.dist_data = copy.copy(dist_data) # history to history distance data used to determine which histories should be added into the alignment next
        self.num_events = 0
        self.num_gaps = 0
    def get_bin(self, index):
        if index > len(self.bin_columns) - 1:
            return False
        return self.bin_columns[index]
    def get_num_gaps(self):
        return self.num_gaps
    def extend_alignment(self, index):
        """
        Extends alignment with empty bin columns until self.bin_columns is as long as index.

        Parameter:
        index (int): how long self.bin_columns should be after extending
        """
        while index > len(self.bin_columns)-1:
            self.bin_columns.append({})
            self.bin_columns[-1]["-"] = self.num_histories
            self.num_gaps += self.num_histories
            self.length += 1

    def add_event_name_to_bin(self, index, event_name, events_in_bin):
        """
        Adds an event_name to self.bin_columns

        Parameters:
        index (int): bin index at which to add event_name
        event_name (str): name of event to add to self.bin_columns
        events_in_bin (int): number of events in original bin
        """
        # TODO: figure out why we divide by the events in a bin here
        # if bin has not been created, create it
        if index > len(self.bin_columns)-1:
            
            self.extend_alignment(index)

        if event_name in self.get_bin(index).keys():
            self.get_bin(index)[event_name] += 1 # TODO: should i remove the divison by events_in_bin?
            
        else:
            self.get_bin(index)[event_name] = 1
        
        self.num_events += 1/events_in_bin

    def add_history(self, history:BinnedHistory, patient:Patient):
        """
        Adds an aligned binned history to self.bin_columns and self.alignments

        Paramters:
        history (BinnedHistory): aligned binned history to add to alignments and self.bin_columns.
        patient (Patient): patient object connected with the history parameter
        """
        # remove patient similarity scores from dist_data
        # TODO: check that this works
        # patient2 is the next_patient so we want to get rid of that option
        

        self.alignments.add_alignment(patient, history)
        self.patients.add(patient)
        for index in sorted(list(history.bins.keys())):
            # bin is empty = indel
            if len(history.bins[index].event_names) == 0:
                self.add_event_name_to_bin(index, "-", 1)
                self.num_gaps += 1
            else:
                for event_name in history.bins[index].event_names:
                    self.add_event_name_to_bin(index, event_name, len(history.bins[index].event_names))
        
        self.num_histories += 1

    def update_dist_data(self):
        """
        Updates self.distance data by removing all data for patients in self.patients from self.dist_data 
        where 'patient 2' == a patient from self.patients. This ensures we don't add duplicate treatment sequences 
        to the profile.
        """
        indices_to_pop = []
        for i, row in enumerate(self.dist_data):
            if row['patient2'] in self.patients:
                indices_to_pop.append(i)
        
        indices_to_pop = reversed(indices_to_pop)
        for i in indices_to_pop:
            self.dist_data.pop(i)
    
    def get_next_patient(self, num_sub_starts):
        """
        Get "num_sub_starts" number of patients that are the most similar to patients 
        already in alignment. 
        """
        self.update_dist_data()
        

        my_data = []    # list of rows where 'patient1' is a patient in self.patients
        for row in self.dist_data:
            if row['patient1'] in self.patients:
                my_data.append(row)

        if num_sub_starts > len(my_data):
            num_sub_starts = len(my_data)

        my_data = sorted(my_data, key=lambda d: d['distance'], reverse=False)
        patients = set()
        for i in range(0, num_sub_starts):
            patients.add(my_data[i]['patient2'])
            if my_data[i]['patient1'] not in self.patients:
                print("error!")
            if my_data[i]['patient2'] in self.patients:
                print("error!")

        return patients
    
    def get_value(self, index, event_name):
        # gets count of an event_name at bin at index
        
        if index > len(self.bin_columns) - 1:
            print("index out of range")
            return False
        
        if event_name not in self.bin_columns[index]:
            return 0 # TODO: update if we want to use psuedocounts
        
        fraction = self.bin_columns[index][event_name] 
        return fraction
    def get_total_events(self, index):
        total_count = 0
        for event_name, count in self.get_bin(index).items():
            if event_name != "-":
                total_count += count
        return total_count
    def get_min_dist(self, MATCH_REWARD, MISMATCH_PENALTY, INDEL_PENALTY):        
        min_dist = 0
        for index in range(self.get_length()):
            total_events = self.get_total_events(index)
            indel_count = self.get_bin(index)["-"]
            indel_dist = indel_count * -INDEL_PENALTY + total_events * MISMATCH_PENALTY
            match_dist = indel_count * INDEL_PENALTY + total_events * MATCH_REWARD
            min_dist += min(indel_dist, match_dist)
        return min_dist
    
    def get_max_dist(self, MATCH_REWARD, MISMATCH_PENALTY, INDEL_PENALTY):
        # things to try: nall events, mismatched event
        max_dist = 0
        for index in range(self.get_length()):
            total_events = self.get_total_events(index)
            indel_count = self.get_bin(index)["-"]
            mismatch_dist = indel_count * INDEL_PENALTY + total_events * MISMATCH_PENALTY
            max_dist += mismatch_dist
        return max_dist
    def get_indel_penalty(self, index, MISMATCH_PENALTY, INDEL_PENALTY):
        penalty = 0
        for event_name, count in self.get_bin(index).items():
            if event_name == "-":
                penalty += count * INDEL_PENALTY
            else:
                penalty += count * MISMATCH_PENALTY
        return penalty


    def get_distance(self, index, bin:Bin, MATCH_REWARD, MISMATCH_PENALTY, INDEL_PENALTY):
        # check that we're not outside the bounds of the profile
        distance = 0
        if index > len(self.bin_columns) - 1:
            return False
        # if bin is empty, we count it as an indel ("-")
        if len(bin.event_names) == 0:
            total_events = self.get_total_events(index)
            indel_count = self.get_value(index, "-")
            distance += (indel_count + 1) * -INDEL_PENALTY 
            distance += total_events * MISMATCH_PENALTY
            
            # fix this later to include mismatch penalty for events missed in profile bin
            # TODO: fix this
            return distance
        
        # add up count_event/num seqs / portion of the bin it takes up
       

        for event_name, count in self.get_bin(index).items():
            
            if event_name in bin.event_names:
                distance += count * MATCH_REWARD
            else:
                if event_name == "-":
                    distance += count * -INDEL_PENALTY
                else:
                    distance += count * MISMATCH_PENALTY
        
        for event_name in bin.event_names:
            if event_name in self.get_bin(index).keys():
                distance += MATCH_REWARD
            else:
                distance += MISMATCH_PENALTY
        
        


            # distance -= self.get_value(index, event_name)/len(bin.event_names)
            
        return distance
    
    def update_alignment(self, profile_indicies):
        self.alignments.update_alignment(profile_indicies)
        new_length = len(profile_indicies)
        if new_length > self.length:
            for i, index in enumerate(profile_indicies):
                if index == "-":
                    self.bin_columns.insert(i, {"-":(self.num_histories)})
                    self.num_gaps += self.num_histories
    
    def get_length(self):
        return len(self.bin_columns)
    
    def print(self):
        for i, bin in enumerate(self.bin_columns):
            s = f"{i}: {bin}"
            print(s)
    def get_simplified_string(self, threshold=0.2):
        s = ""
        for i, bin in enumerate(self.bin_columns):
           
                
            # get total num events and indel
            total_count = 0
            for event_name, count in bin.items():
                total_count += count
            
            
            
            s_ = "["
            j = 0
            for event_name, count in bin.items():
                
                # if s_ != "[":
                #     s_ += ", "
                if count/total_count >= threshold:
                    if event_name != "-":
                        if j != 0:
                            s_ += ", "
                        if count/total_count != 1.0:
                            s_ += f"{event_name}: {round(count/total_count, 2)}"
                        else:
                            s_ += f"{event_name}"
                        j += 1
                
            s_ += "]"
            if len(s_) != 2:
                s += s_
                s += " "
        return s
    
    def get_simplified_dict(self, threshold=0.05, med_letters=3):
        data = {}

        # figure out length of bins
        total_bins = 0
        bins_to_keep = []
        for i, bin in enumerate(self.bin_columns):
                           
            # get total num events and indel
            total_count = 0

            for event_name, count in bin.items():
                total_count += count

            for event_name, count in bin.items():

                if count/total_count >= threshold:

                    if event_name != "-":
                        total_bins += 1
                        bins_to_keep.append(bin)
                        break
 
            
        for i, bin in enumerate(bins_to_keep):
                           
            # get total num events and indel
            total_count = 0
            for event_name, count in bin.items():
                total_count += count

            # add frequencies to data
            for event_name, count in bin.items():
                
                if count/total_count >= threshold:
                    if event_name != "-":
                        # add event to data
                        abbrv = self.get_med_abbrv(event_name, med_letters)
                        if abbrv not in data.keys():
                            data[abbrv] = [0]*total_bins

                        data[abbrv][i] = count/(self.num_histories)
        return data

    def get_med_abbrv(self, med_name, med_letters):
        if len(med_name) >= med_letters:
            abbrv = med_name[:med_letters].upper()
            
        else:
            abbrv = copy.copy(med_name)
            while len(abbrv) < med_letters:
                abbrv += "-"
        
        return abbrv
        
    def get_logo(self, file_path, color_dict=None, threshold=0.2, title=None, med_letters=3):
        data = self.get_simplified_dict(threshold, med_letters)
        print(len(data))
        if len(data) == 0:
            print("consensus sequence length = 0")
            return
        if color_dict is None:
            logo = Logo(pd.DataFrame(data), baseline_width=0.2)
        else:
            logo = Logo(pd.DataFrame(data), color_scheme=color_dict, baseline_width=0.2)
        # Style the logo
        logo.style_spines(visible=False)
        logo.style_spines(spines=['left', 'bottom'], visible=True)
        logo.style_xticks(rotation=90, fmt='%d', anchor=0)
        logo.ax.set_ylabel('Frequency')
        # logo.ax.set_ylim([0, 1])
        if title:
            logo.ax.set_title(title)
        logo.fig.savefig(file_path)
        


                
                

    
