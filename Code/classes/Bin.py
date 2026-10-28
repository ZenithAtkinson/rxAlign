from Event  import *

class Bin(object):
    """
    Data structure for storing events that occur within the bin's timeframe.
    """
    def __init__(self):
        """
        Construct a new 'Bin' object.
        """
        self.event_names = set()      # Set of names of events. This will need to be changed to be a list of 
                                      # event objects if we want to store more information about the event (event type, med_class, dosage, etc)

    def add_event(self, event:Event):
        """
        Adds an event.
        """
        self.event_names.add(event.name)

    def get_event(self, event_name):  # TODO: update this later to return event itself
        """
        Checks if event exists in bin.
        
        Returns:
        event_name if it exists, otherwise returns False
        """
        if event_name in self.event_names:
            return event_name
        else: 
            return False
        
    def get_distance(self, bin, MATCH_REWARD, MISMATCH_PENALTY, INDEL_PENALTY):
        """
        Compares self with another bin.

        Parameters:
        bin (Bin): Bin to compare self to
        MATCH_REWARD (int or float): Reward for matching events. Should be negative.
        MISMATCH_PENALTY (int or float): Penalty for an event not occuring in both bins. Should be positive.
        INDEL_PENALTY (int or float): -INDEL_PENALTY is used as the distance of two empty bins. Should be positive.
        
        Returns:
        Distance between the two bins
        """

        distance = 0
        
        if len(self.event_names) == 0:
            if len(bin.event_names) == 0:
                return -INDEL_PENALTY       # this could be changed to allow more user control
        

            
        # TODO: update this to work compare category, not just names of events.
        
        
        total_events = len(self.event_names) + len(bin.event_names)
        unmatched_events = total_events 
        
        for event_name in bin.event_names:
            if event_name in self.event_names:
                unmatched_events = unmatched_events - 2
        
        distance = unmatched_events*MISMATCH_PENALTY + (total_events - unmatched_events) * MATCH_REWARD
        return distance
    
    def to_string(self):
        """
        Converts bin into string containing the bin's events.
        """
        s = "["
        for i, event in enumerate(self.event_names):
            s += event
            if i != len(self.event_names)-1:
                s +=  ", "
        s += "]"
        return s
