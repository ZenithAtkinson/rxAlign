from Timeline import Timeline

#jaccard algorithm: for comparing two sets (one event vs one event) - event_cost
#needleman-wunsch: two sequences (a whole timeline vs a whole timeline) - distance

#1: Calculate how different two events are, on a scale of 0 to 1
def event_cost(e1, e2):
    e1 = set(e1.treatments.treatment)
    e2 = set(e2.treatments.treatment)
    shared = len(e1 & e2)
    total = len(e1 | e2)
    return (1 - (shared / total))
#2: Get the dynamic programming table
def distance(t1, t2, indel=1.0):
    seq1 = t1.event_list
    seq2 = t2.event_list
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    val_matrix = [[0] * (len_seq2+1) for _ in range(len_seq1+1)]

    #init the first row/col
    for i in range(len_seq1 + 1):
        val_matrix[i][0] = i * indel
    for j in range(len_seq2 + 1):
        val_matrix[0][j] = j * indel
    
    #fill the rest of the table: each cell = cheapest of match, delete, insert
    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            match = val_matrix[i-1][j-1] + event_cost(seq1[i-1], seq2[j-1])
            delete = val_matrix[i-1][j] + indel   # skip an event from seq1
            insert = val_matrix[i][j-1] + indel   # skip an event from seq2
            val_matrix[i][j] = min(match, delete, insert)

    #bottom-right cell = raw distance between the whole timelines
    return val_matrix[len_seq1][len_seq2]

#3: Normalize! Scales the raw distance to 0 (identical) .. 1 (as different as possible for these lengths)
def normalized_distance(t1, t2, indel=1.0):
    len_seq1 = len(t1.event_list)
    len_seq2 = len(t2.event_list)

    # worst case: each paired-up event is either a full mismatch (cost 1) or a
    # delete + insert (cost 2*indel), whichever is cheaper; leftover events are indels.
    # with indel = 1 this is just max(len_seq1, len_seq2)
    paired = min(len_seq1, len_seq2)
    leftover = abs(len_seq1 - len_seq2)
    max_distance = paired * min(1, 2 * indel) + leftover * indel

    if max_distance == 0:   # both timelines empty -> identical
        return 0.0

    return distance(t1, t2, indel) / max_distance


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    from Event import Event
    from Treatment import Treatment

    # helper: make("A", "B") -> timeline A -> B ; use a list for combos: make(["A", "B"], "C") -> A+B -> C
    def make(*events):
        t = Timeline("test")
        for start, drugs in enumerate(events):
            if isinstance(drugs, str):
                drugs = [drugs]
            t.add_event(Event(Treatment(list(drugs)), start, start + 1))
        return t

    # --- event_cost examples ---
    ev_AB = Event(Treatment(["A", "B"]), 0, 1)
    ev_A  = Event(Treatment(["A"]), 0, 1)
    ev_BC = Event(Treatment(["B", "C"]), 0, 1)
    ev_C  = Event(Treatment(["C"]), 0, 1)

    print("event_cost")
    print("  {A,B} vs {A,B}:", event_cost(ev_AB, ev_AB), "  expect 0")
    print("  {A,B} vs {A}:  ", event_cost(ev_AB, ev_A),  "  expect 0.5")
    print("  {A,B} vs {B,C}:", round(event_cost(ev_AB, ev_BC), 3), "  expect 0.667")
    print("  {A} vs {C}:    ", event_cost(ev_A, ev_C),   "  expect 1")

    # --- timeline examples ---
    t_AB      = make("A", "B")             # A -> B
    t_A       = make("A")                  # A
    t_B       = make("B")                  # B
    t_ABC     = make("A", "B", "C")        # A -> B -> C
    t_AC      = make("A", "C")             # A -> C
    t_ABcombo = make(["A", "B"], "C")      # A+B -> C
    t_BA      = make("B", "A")             # B -> A
    t_empty   = Timeline("empty")          # no events

    # (name, t1, t2, expected raw, expected normalized)
    tests = [
        ("A->B   vs A->B  ", t_AB,      t_AB,    0,   0),
        ("A->B   vs A     ", t_AB,      t_A,     1,   0.5),
        ("A      vs B     ", t_A,       t_B,     1,   1),
        ("A->B->C vs A->C ", t_ABC,     t_AC,    1,   0.333),
        ("A+B->C vs A->C  ", t_ABcombo, t_AC,    0.5, 0.25),
        ("A->B   vs B->A  ", t_AB,      t_BA,    2,   1),
        ("empty  vs A     ", t_empty,   t_A,     1,   1),
        ("empty  vs empty ", t_empty,   t_empty, 0,   0),
    ]

    print("\ndistance (raw / normalized)")
    for name, t1, t2, exp_raw, exp_norm in tests:
        raw = distance(t1, t2)
        norm = round(normalized_distance(t1, t2), 3)
        symmetric = distance(t1, t2) == distance(t2, t1)
        ok = "OK " if (raw == exp_raw and norm == exp_norm and symmetric) else "BAD"
        print(f"  {ok} {name} raw={raw} (expect {exp_raw})  norm={norm} (expect {exp_norm})  symmetric={symmetric}")
