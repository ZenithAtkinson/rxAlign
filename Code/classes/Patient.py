from BinnedHistory import *
from Event import *

class Patient(object):
    """
    Data structure for storing a patient's id and binned history.
    """
    def __init__(self, id, bin_length):
        """
        Construct a new patient object and initialize patient's binned history.
        """
        self.id = id
        # self.age_at_diagnosis = age_at_diagnosis
        self.binned_history = BinnedHistory(bin_length)
        
    def add_event(self, event:Event):
        """
        Add event to patient's binned history
        """
        self.binned_history.add_event(event)
