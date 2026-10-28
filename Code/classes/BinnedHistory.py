from Event import *
from Bin  import *

class BinnedHistory(object):
    """
    Stores a treatment path in the form of an index:Bin dictionary 
    """
    def __init__(self, bin_length=0, bins:list=None):
        """
        Construct a new 'BinnedHistory' object. Will populate 'BinnedHistory' using bins parameter if provided. 

        Parameters:
        bin_length (int): length (in time) of all bins
        bins (list): chronological list of 'Bins' objects
        """
        self.bin_length = bin_length
        self.bins = {}                  # key=index, value=Bin
        self.last_bin_index = 0
        
        if bins:
            self.create_binned_history_from_list(bins) 

    def get_num_events(self):
        """
        Returns total number of events.
        """
        num_events = 0
        for bin in self.bins.values():
            num_events += len(bin.event_names)
        
        return num_events
            
    def create_binned_history_from_list(self, bins):
        """
        Iteratively inserts bins into self.bins

        Parameter:
        bins (list): chronological list of 'Bin' objects
        """
        self.last_bin_index = -1
        for bin in bins:
            self.append_bin(bin)
            
    def to_string(self):
        """
        Returns a string representation of self.bins
        """
        s = ""
        for index in sorted(list(self.bins.keys())):
            s += f'{self.bins[index].to_string()}'
        return s
    
    def update_last_bin_index(self, new_index):
        """
        Updates self.last_bin_index to new_index
        """
        self.last_bin_index = new_index

    def get_length(self):
        """
        Returns length of 'BinnedHistory'
        """
        return self.last_bin_index + 1
    
    def add_event(self, event:Event):
        """
        Adds an event to self.bins. Note: event will be added to multiple bins
        if its duration is longer than self.bin_length

        Parameter:
        event (Event): event to add to binned history
        """
        start_index = int(event.start // self.bin_length)
        end_index = int(event.end // self.bin_length)
        
        index = start_index
        while index <= end_index:
            bin = self.get_bin(index)

            # if bin hasn't been created yet, make a new one
            if bin == False:
                bin = self.add_bin(index)

            bin.add_event(event)
            index += 1

    def add_bin(self, index):
        """
        Create a new bin at the given index in self.bins.
        """
        bin = Bin()
        self.bins[index] = bin
        if index > self.last_bin_index:
            self.update_last_bin_index(index)
        return bin
    
    def get_bin(self, index):
        """
        Returns the bin at the specified index. If index is outside of binned history range, returns False
        """
        if type(index) == str:
            index = int(index) # after loading objects with pickle indices that used to be ints are now strings
        
        if index in self.bins.keys():
            return self.bins[index]
        
        else:
            # bin exists within the bounds of the patient timeline
            if index < self.last_bin_index:
                bin = self.add_bin(index)
                return bin
            
            else:
                return False
            
    def append_bin(self, bin):
        """
        Appends a bin to the end of the binned history.
        """       
        self.bins[self.last_bin_index+1] = bin
        self.last_bin_index += 1

    def get_num_gaps(self):
        """
        Returns the number of bins without events in the binned history.
        """
        num_gaps = 0
        index = 0

        while index <= self.last_bin_index:

            bin = self.get_bin(index)

            if bin == False:
                # gap
                num_gaps += 1

            elif len(bin.event_names) == 0:
                num_gaps += 1
           
            index += 1
            
        return num_gaps
