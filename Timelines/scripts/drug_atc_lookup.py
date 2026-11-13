# Looks up the WHO ATC class of every drug in the GENIE regimen dataset using NLM's RxClass API,
# turns each ATC code into a MEC-TX treatment class, and compares that to the current mapping
# (Code/data/med_class.tsv + the overrides in Timelines/R_scripts/mectx_run.R).
#
# Run:     python Timelines/scripts/drug_atc_lookup.py      (from anywhere; needs internet)
# Output:  Timelines/scripts/drug_atc_classes.csv           (one row per drug)
#
# Only drug names are sent to the API (https://rxnav.nlm.nih.gov), never patient data.
# API docs: https://lhncbc.nlm.nih.gov/RxNav/APIs/RxClassAPIs.html
# ATC index (to look up any code by hand): https://atcddd.fhi.no/atc_ddd_index/

import csv
import json
import os
import time
import urllib.parse
import urllib.request
from collections import Counter

# ---------- paths (found relative to this file, so it works on the laptop and at home) ----------
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA_PATH = os.path.join(ROOT, "Data", "regimen_cancer_level_dataset.csv")
MED_CLASS_PATH = os.path.join(ROOT, "Code", "data", "med_class.tsv")
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "drug_atc_classes.csv")

API = "https://rxnav.nlm.nih.gov/REST/rxclass/class/byDrugName.json"

# ---------- ATC code -> MEC-TX class ----------
# Checked in order; the first matching prefix wins, so specific codes go before general ones.
# "REVIEW" = ATC alone can't decide; look at the drug by hand.
ATC_RULES = [
    ("L01FF", "IO"),              # PD-1 / PD-L1 inhibitors (nivolumab, pembrolizumab, atezolizumab, ...)
    ("L01A",  "Chemo"),           # alkylating agents
    ("L01B",  "Chemo"),           # antimetabolites
    ("L01C",  "Chemo"),           # plant alkaloids (taxanes, vinca alkaloids, etoposide)
    ("L01D",  "Chemo"),           # cytotoxic antibiotics
    ("L01XA", "Chemo"),           # platinum compounds
    ("L01XB", "Chemo"),           # methylhydrazines (procarbazine)
    ("L01E",  "Small_Molecule"),  # protein kinase inhibitors (EGFR, ALK, BRAF, MEK, CDK4/6, mTOR, ...)
    ("L01XG", "Small_Molecule"),  # proteasome inhibitors
    ("L01XH", "Small_Molecule"),  # histone deacetylase inhibitors
    ("L01XJ", "Small_Molecule"),  # hedgehog pathway inhibitors
    ("L01XK", "Small_Molecule"),  # PARP inhibitors
    ("L01F",  "Targeted"),        # monoclonal antibodies + antibody-drug conjugates (bevacizumab, trastuzumab, ...)
    ("L01XX", "REVIEW"),          # "other antineoplastic agents": a mixed bag (arsenic trioxide, hydroxyurea, ...)
    ("L01",   "REVIEW"),          # anything else in L01
    ("L02",   "Hormone"),         # endocrine therapy
    ("H01CB", "Hormone"),         # somatostatin analogues (lanreotide, octreotide; used for neuroendocrine tumours)
    ("L03AA", "Ancillary"),       # colony stimulating factors (filgrastim): supportive care
    ("L03",   "IO"),              # immunostimulants: interferons, interleukins, BCG
    ("L04",   "REVIEW"),          # immunosuppressants (includes lenalidomide / pomalidomide)
    ("V03AF", "Ancillary"),       # detoxifying agents for antineoplastic treatment (leucovorin, mesna)
    ("H02",   "Ancillary"),       # corticosteroids
    ("A04",   "Ancillary"),       # antiemetics
    ("M05",   "Ancillary"),       # bone drugs (zoledronic acid, denosumab)
]

# When a drug has several ATC codes (e.g. megestrol: G03AC and L02AB), prefer the cancer-relevant one.
ATC_PRIORITY = ["L01", "L02", "L03", "L04", "V03", "H02", "A04", "M05"]

# Drugs whose ATC code gives the wrong MEC-TX class.
NAME_EXCEPTIONS = {
    "IPILIMUMAB":   ("IO", "CTLA-4 checkpoint inhibitor; ATC files it under L01FX (other antibodies)"),
    "TREMELIMUMAB": ("IO", "CTLA-4 checkpoint inhibitor; ATC files it under L01FX (other antibodies)"),
}

# Drugs RxClass has no ATC code for (not in RxNorm, or a placeholder name). These are judgement calls,
# NOT from the API: verify each one (NCI Drug Dictionary: https://www.cancer.gov/publications/dictionaries/cancer-drug)
# and get your PI to sign off. They show up with status "MANUAL" in the output.
MANUAL = {
    "ENSARTINIB":           ("Small_Molecule", "ALK tyrosine kinase inhibitor; not FDA-approved during data collection, so not in RxNorm"),
    "ICOTINIB":             ("Small_Molecule", "EGFR tyrosine kinase inhibitor; approved in China only, so not in RxNorm"),
    "ROCILETINIB":          ("Small_Molecule", "EGFR tyrosine kinase inhibitor; development stopped, not in RxNorm"),
    "NAVITOCLAX":           ("Small_Molecule", "BCL-2/BCL-xL inhibitor; investigational, not in RxNorm"),
    "BCG SOLUTION":         ("IO", "WHO ATC lists BCG vaccine as L03AX03 (immunostimulant), but RxClass has no ATC link for it"),
    "BCG VACCINE":          ("IO", "WHO ATC lists BCG vaccine as L03AX03 (immunostimulant), but RxClass has no ATC link for it"),
    "INVESTIGATIONAL DRUG": ("Others", "placeholder name: the actual drug is unknown"),
    "OTHER ANTINEOPLASTIC": ("Others", "placeholder name: the actual drug is unknown"),
    "OTHER NOS":            ("Others", "placeholder name: the actual drug is unknown"),
    "OTHER HORMONE":        ("Hormone", "placeholder name, but explicitly a hormone therapy"),
}

# Names RxClass doesn't recognise -> the name to query instead (formulations of a known drug).
ALIASES = {
    "NABPACLITAXEL": "paclitaxel",
    "PACLITAXEL LOADED POLYMERIC MICELLE": "paclitaxel",
    "PACLITAXEL POLIGLUMEX": "paclitaxel",
    "PACLITAXEL TREVATIDE": "paclitaxel",
    "PEGYLATED LIPOSOMAL DOXORUBICIN HYDROCHLORIDE": "doxorubicin",
    "VINCRISTINE SULFATE LIPOSOME": "vincristine",
    "RITUXIMAB AND HYALURONIDASE HUMAN": "rituximab",
    "TRASTUZUMAB EMTANSINE": "ado-trastuzumab emtansine",
    "BCG SOLUTION": "BCG, live",
    "BCG VACCINE": "BCG, live",
}

# Salt / ester words stripped off as a last attempt ("CABOZANTINIB SMALATE" -> "cabozantinib").
SALT_WORDS = {"HYDROCHLORIDE", "DISODIUM", "SODIUM", "ACETATE", "MESYLATE", "TARTRATE", "SULFATE",
              "CITRATE", "PHOSPHATE", "CALCIUM", "DIMALEATE", "SMALATE", "MALATE", "TOSYLATE"}

# ---------- current mapping: copy of what mectx_run.R does ----------
TSV_RENAME = {"Targeted": "Targeted", "Chemo": "Chemo", "ICB": "IO", "Other": "Others"}
OVERRIDES = {   # keep in sync with the overrides list in Timelines/R_scripts/mectx_run.R
    "Targeted": ["AFATINIB DIMALEATE", "CABOZANTINIB SMALATE", "DARATUMUMAB", "ENSARTINIB", "ENTRECTINIB",
                 "ERLOTINIB HYDROCHLORIDE", "ICOTINIB", "NAVITOCLAX", "NECITUMUMAB", "RITUXIMAB AND HYALURONIDASE HUMAN",
                 "ROCILETINIB", "RUXOLITINIB PHOSPHATE", "TRASTUZUMAB EMTANSINE"],
    "Chemo": ["ARSENIC TRIOXIDE", "BLEOMYCIN", "BUSULFAN", "CHLORAMBUCIL", "ETOPOSIDE PHOSPHATE", "LURBINECTEDIN",
              "MECHLORETHAMINE HYDROCHLORIDE", "PACLITAXEL LOADED POLYMERIC MICELLE", "PACLITAXEL POLIGLUMEX",
              "PACLITAXEL TREVATIDE", "PEGYLATED LIPOSOMAL DOXORUBICIN HYDROCHLORIDE", "VINCRISTINE SULFATE LIPOSOME"],
    "IO": ["TREMELIMUMAB", "BCG SOLUTION", "BCG VACCINE"],
    "Hormone": ["ANASTROZOLE", "LETROZOLE", "TAMOXIFEN CITRATE", "EXEMESTANE", "LEUPROLIDE ACETATE", "BICALUTAMIDE",
                "DEGARELIX", "MEGESTROL ACETATE", "ABIRATERONE ACETATE", "LANREOTIDE ACETATE", "OTHER HORMONE"],
    "Ancillary": ["LEUCOVORIN CALCIUM"],
    "Others": ["INVESTIGATIONAL DRUG", "OTHER ANTINEOPLASTIC", "OTHER NOS"],
}


def current_mapping():
    mapping = {}
    with open(MED_CLASS_PATH, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            mapping[row["Medication"].strip().upper()] = TSV_RENAME[row["Class"].strip()]
    for cls, names in OVERRIDES.items():
        for name in names:
            mapping[name] = cls
    return mapping


def drugs_in_data():
    """Base drug name -> (records in all cancers, records in index cancers), same filters as mectx_run.R."""
    all_n, index_n = Counter(), Counter()
    with open(DATA_PATH, newline="") as f:
        for row in csv.DictReader(f):
            for n in range(1, 6):
                drug = row[f"drugs_drug_{n}"]
                if drug == "" or row[f"dx_drug_start_int_{n}"] == "":
                    continue
                base = drug.split("(", 1)[0].strip().upper()   # drop the "(synonyms...)"
                all_n[base] += 1
                if row["redcap_ca_index"] == "Yes":
                    index_n[base] += 1
    return all_n, index_n


def query_atc(name):
    """Returns the sorted list of (ATC code, ATC class name) RxClass gives for a drug name."""
    url = API + "?" + urllib.parse.urlencode({"drugName": name, "relaSource": "ATC"})
    with urllib.request.urlopen(url, timeout=30) as r:
        data = json.load(r)
    time.sleep(0.1)   # RxNav asks for at most 20 requests / second
    items = data.get("rxclassDrugInfoList", {}).get("rxclassDrugInfo", [])
    # RxClass also returns combination products that contain the drug (tty "MIN", e.g. abiraterone/niraparib
    # -> L01XK PARP inhibitors). Ignore those unless the drug only exists as a combination.
    single = [i for i in items if i["minConcept"]["tty"] != "MIN"]
    items = single or items
    return sorted({(i["rxclassMinConceptItem"]["classId"], i["rxclassMinConceptItem"]["className"]) for i in items})


def lookup(drug):
    """Tries the name as-is, then an alias, then the name without salt words. Returns (query used, codes)."""
    attempts = [drug.lower()]
    if drug in ALIASES:
        attempts.append(ALIASES[drug])
    stripped = " ".join(w for w in drug.split() if w not in SALT_WORDS)
    if stripped and stripped != drug:
        attempts.append(stripped.lower())
    for q in attempts:
        codes = query_atc(q)
        if codes:
            return q, codes
    return attempts[0], []


def atc_to_class(codes):
    """Picks the most cancer-relevant ATC code and converts it with ATC_RULES."""
    if not codes:
        return "", "REVIEW"
    def priority(code):
        return next((i for i, p in enumerate(ATC_PRIORITY) if code.startswith(p)), len(ATC_PRIORITY))
    best = min((c for c, _ in codes), key=priority)
    cls = next((c for prefix, c in ATC_RULES if best.startswith(prefix)), "REVIEW")
    return best, cls


def main():
    current = current_mapping()
    all_n, index_n = drugs_in_data()
    rows = []
    for i, drug in enumerate(sorted(all_n), 1):
        query, codes = lookup(drug)
        best, atc_class = atc_to_class(codes)
        note = ""
        if drug in NAME_EXCEPTIONS:
            atc_class, note = NAME_EXCEPTIONS[drug]
        cur = current.get(drug, "UNMAPPED")

        if not codes and drug in MANUAL:
            atc_class, note = MANUAL[drug]
            status = "MANUAL: no ATC code, judgement call (verify)"
        elif not codes:
            status = "REVIEW: no ATC code found"
        elif atc_class == "REVIEW":
            status = f"REVIEW: ATC {best} doesn't decide"
        elif atc_class == cur:
            status = "agree"
        elif atc_class == "Small_Molecule" and cur == "Targeted":
            status = "agree (Targeted -> Small_Molecule split)"
        else:
            status = "DISAGREE"

        rows.append({
            "drug": drug,
            "records_all": all_n[drug],
            "records_index": index_n[drug],
            "current_class": cur,
            "atc_class": atc_class,
            "status": status,
            "atc_code_used": best,
            "all_atc_codes": "; ".join(f"{c} {n}" for c, n in codes),
            "query_used": query,
            "note": note,
        })
        print(f"[{i}/{len(all_n)}] {drug}: {best or '-'} -> {atc_class} ({status})")

    with open(OUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} drugs to {OUT_PATH}")
    for status, n in Counter(r["status"].split(":")[0] for r in rows).most_common():
        print(f"  {status}: {n}")


if __name__ == "__main__":
    main()
