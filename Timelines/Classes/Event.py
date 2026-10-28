# A point in time on the timeline where one or mulitple treatments were administered.

# Will contain a treatment list (sorted alphabetically), start date, end date, index in the timeline, 
# and what other events it points to (ahead/AND behind?)
from Treatment import *

class Event:
    def __init__(self, treatments, start, end, ):
        self.treatments = treatments
        self.start = start
        self.end = end

me = Treatment(["organe", "lemon"])

lesgo = Event(me, 1, 2)

print(lesgo.start)
