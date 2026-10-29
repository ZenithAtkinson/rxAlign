# A timeline for a specified patient. Can be stored in a list of timelines of all other patients.
# The timeline will be a sequence of events, each one pointed to the next. (each event has treatments and start/end date for them)
# Each timeline will also have an identifier that is the same as the patient identifier

# the events in the timeline will be sorted according to their *START* date. Needs a function to add, remove, etc events to a given timeline
# Data is being pulled from ~/data/regimen_cancer_level_dataset.csv

# Required variables:
# timeline object, which is a list of events (the events are sorted by start date, and each event is then also assigned an index(?) and a "points_to" tag(?))
#? = are these even needed or just unnecessary?
# patient ID
# 

# Required functions:
# Add event (given a start date, which it uses to place it in the timeline and point it to corresponding events)
    # will we ever have to remove an event? I dont think so
# Assign a patient ID

class Timeline:
    def __init__(self, patient_ID=None):
        self.patient_ID = patient_ID if patient_ID is not None else "Unassigned"
        self.event_list = []
        pass

    




