import csv
from Timeline import Timeline 
from Event import Event
from Treatment import Treatment

def load_timelines(path):
    timelines = {}
    with open(path, newline="") as f:
        # looking at datarows, if patient doesnt exist, add it to the timeline dictn
        for row in csv.DictReader(f):
            patient = row["record_id"]
            if patient not in timelines:
                timelines[patient] = Timeline(patient)
            patient_data = timelines[patient]

            for n in range(1, 6):              
                drug = row[f"drugs_drug_{n}"]
                if drug == "":                 # this regimen has fewer than n drugs
                    continue                   # skip to the next n
                # clean the drug name, read start/end, make an Event, add_event it
                start = row[f"drugs_startdt_int_{n}"]
                if start == "":
                    continue
                start = float(start) #making it work for any "date" columns we want to use

                end = row[f"dx_drug_end_int_{n}"]
                if end == "":
                    continue
                end = float(end) #making it work for any "date" columns we want to use

                #adding these to a new event:
                patient_data.add_event(Event(Treatment([drug]), start, end))


    return timelines

timelines = load_timelines("C:/Users/Zenith/PhDCoding/rxAlign/Data/regimen_cancer_level_dataset.csv")
# Testing first patient in timeline
first = list(timelines.values())[0]
second = list(timelines.values())[1]
third = list(timelines.values())[2]
print(first.patient_ID)
for e in first.event_list:
    print(e.start, e.treatments)

print(second.patient_ID)
for e in second.event_list:
    print(e.start, e.treatments)

print(third.patient_ID)
for e in third.event_list:
    print(e.start, e.treatments)
print(third.event_list)