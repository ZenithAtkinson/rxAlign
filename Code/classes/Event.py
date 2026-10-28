class Event(object):
    """
    Stores an event and its information
    """
    def __init__(self, name, start, end, category=False):
        """
        Constructs an 'Event' object.

        Parameters:
        name (str): event name (in the case of a medication, this is the medication name)
        start (float or int): event start time
        end (float or int): event end time
        category: event category
        """
        self.name = name
        self.start = start
        self.end = end
        self.category = category # TODO: implement this once we want to use med class or event type