# Given to the event class, this will be the list of treatments sorted alphabetically. Putting it into its own class like this isnt
# totally necessary, but it gives more flexibility and allows for future additions more easily to the overall timeline

class Treatment:
    def __init__(self, treatment=None):
        self.treatment = treatment if treatment is not None else []
        self.length = len(treatment) if treatment is not None else 0

    def __repr__(self):
    # Just did this so if you print "treatment" by itself, itll return the list of treatments
        return str(self.treatment)

    def __str__(self):
        return str(self.treatment)