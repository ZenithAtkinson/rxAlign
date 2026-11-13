# Matches every drug a patient has to a treatment category, with the reason for that category.
# Drug name = the text before "(", e.g. "Carboplatin(Paraplatin,...)" -> "CARBOPLATIN".
# Each group below = one category + one reason + the drugs it applies to.
# Reasons are the drug's WHO ATC code (looked up by Timelines/scripts/drug_atc_lookup.py via NLM RxClass).
# "JUDGEMENT CALL" = the ATC code doesn't settle the category (or there is none); check these with your PI.
# Stops with an error if any drug in the data isn't in a group below.
#
# Run:     Rscript Timelines/R_scripts/drug_categories.R
# Output:  Timelines/R_scripts/patient_drug_categories.csv  (one row per patient drug: record_id, drug, category, reason)
#
# To undo the Small_Molecule / Targeted split: change category = "Small_Molecule" to "Targeted".

# ALEX - switch these when im remote vs at home
#ROOT <- "C:/Users/default.DESKTOP-UOBTTU9/Bioinformatics Lesgo/rxAlign"
ROOT <- "C:/Users/Zenith/PhDCoding/rxAlign"

groups <- list(
  list(category = "Chemo", reason = "ATC L01AA (Nitrogen mustard analogues): cytotoxic chemotherapy",
       drugs = c("BENDAMUSTINE HYDROCHLORIDE", "CHLORAMBUCIL", "CYCLOPHOSPHAMIDE", "MECHLORETHAMINE HYDROCHLORIDE")),
  list(category = "Chemo", reason = "ATC L01AB (Alkyl sulfonates): cytotoxic chemotherapy",
       drugs = c("BUSULFAN")),
  list(category = "Chemo", reason = "ATC L01AC (Ethylene imines): cytotoxic chemotherapy",
       drugs = c("THIOTEPA")),
  list(category = "Chemo", reason = "ATC L01AD (Nitrosoureas): cytotoxic chemotherapy",
       drugs = c("LOMUSTINE")),
  list(category = "Chemo", reason = "ATC L01AX (Other alkylating agents): cytotoxic chemotherapy",
       drugs = c("DACARBAZINE", "TEMOZOLOMIDE")),
  list(category = "Chemo", reason = "ATC L01BA (Folic acid analogues): cytotoxic chemotherapy",
       drugs = c("METHOTREXATE", "PEMETREXED DISODIUM")),
  list(category = "Chemo", reason = "ATC L01BB (Purine analogues): cytotoxic chemotherapy",
       drugs = c("FLUDARABINE PHOSPHATE")),
  list(category = "Chemo", reason = "ATC L01BC (Pyrimidine analogues): cytotoxic chemotherapy",
       drugs = c("CAPECITABINE", "CYTARABINE", "DECITABINE", "FLUOROURACIL", "GEMCITABINE HYDROCHLORIDE")),
  list(category = "Chemo", reason = "ATC L01CA (Vinca alkaloids and analogues): cytotoxic chemotherapy",
       drugs = c("VINBLASTINE SULFATE", "VINCRISTINE SULFATE", "VINCRISTINE SULFATE LIPOSOME", "VINORELBINE TARTRATE")),
  list(category = "Chemo", reason = "ATC L01CB (Podophyllotoxin derivatives): cytotoxic chemotherapy",
       drugs = c("ETOPOSIDE", "ETOPOSIDE PHOSPHATE")),
  list(category = "Chemo", reason = "ATC L01CD (Taxanes): cytotoxic chemotherapy",
       drugs = c("DOCETAXEL", "NABPACLITAXEL", "PACLITAXEL", "PACLITAXEL LOADED POLYMERIC MICELLE", "PACLITAXEL POLIGLUMEX", "PACLITAXEL TREVATIDE")),
  list(category = "Chemo", reason = "ATC L01CE (Topoisomerase 1 (TOP1) inhibitors): cytotoxic chemotherapy",
       drugs = c("IRINOTECAN HYDROCHLORIDE")),
  list(category = "Chemo", reason = "ATC L01DB (Anthracyclines and related substances): cytotoxic chemotherapy",
       drugs = c("DOXORUBICIN HYDROCHLORIDE", "EPIRUBICIN HYDROCHLORIDE", "PEGYLATED LIPOSOMAL DOXORUBICIN HYDROCHLORIDE")),
  list(category = "Chemo", reason = "ATC L01DC (Other cytotoxic antibiotics): cytotoxic chemotherapy",
       drugs = c("BLEOMYCIN", "MITOMYCIN")),
  list(category = "Chemo", reason = "ATC L01XA (Platinum compounds): cytotoxic chemotherapy",
       drugs = c("CARBOPLATIN", "CISPLATIN", "OXALIPLATIN")),
  list(category = "Chemo", reason = "ATC L01XB (Methylhydrazines): cytotoxic chemotherapy",
       drugs = c("PROCARBAZINE HYDROCHLORIDE")),
  list(category = "Chemo", reason = "JUDGEMENT CALL: ATC L01XX (Other antineoplastic agents) is a catch-all; classed as Chemo because these act as cytotoxics",
       drugs = c("ARSENIC TRIOXIDE", "HYDROXYUREA", "LURBINECTEDIN")),
  list(category = "Small_Molecule", reason = "ATC L01EA (BCR-ABL tyrosine kinase inhibitors): targeted small-molecule inhibitor",
       drugs = c("DASATINIB", "IMATINIB MESYLATE")),
  list(category = "Small_Molecule", reason = "ATC L01EB (Epidermal growth factor receptor (EGFR) tyrosine kinase inhibitors): targeted small-molecule inhibitor",
       drugs = c("AFATINIB DIMALEATE", "ERLOTINIB HYDROCHLORIDE", "GEFITINIB", "OSIMERTINIB")),
  list(category = "Small_Molecule", reason = "ATC L01EC (B-Raf serine-threonine kinase (BRAF) inhibitors): targeted small-molecule inhibitor",
       drugs = c("DABRAFENIB", "VEMURAFENIB")),
  list(category = "Small_Molecule", reason = "ATC L01ED (Anaplastic lymphoma kinase (ALK) inhibitors): targeted small-molecule inhibitor",
       drugs = c("ALECTINIB", "BRIGATINIB", "CERITINIB", "CRIZOTINIB", "LORLATINIB")),
  list(category = "Small_Molecule", reason = "ATC L01EE (Mitogen-activated protein kinase (MEK) inhibitors): targeted small-molecule inhibitor",
       drugs = c("COBIMETINIB", "TRAMETINIB")),
  list(category = "Small_Molecule", reason = "ATC L01EF (Cyclin-dependent kinase (CDK) inhibitors): targeted small-molecule inhibitor",
       drugs = c("PALBOCICLIB")),
  list(category = "Small_Molecule", reason = "ATC L01EG (Mammalian target of rapamycin (mTOR) kinase inhibitors): targeted small-molecule inhibitor",
       drugs = c("EVEROLIMUS")),
  list(category = "Small_Molecule", reason = "ATC L01EH (Human epidermal growth factor receptor 2 (HER2) tyrosine kinase inhibitors): targeted small-molecule inhibitor",
       drugs = c("NERATINIB")),
  list(category = "Small_Molecule", reason = "ATC L01EJ (Janus-associated kinase (JAK) inhibitors): targeted small-molecule inhibitor",
       drugs = c("RUXOLITINIB PHOSPHATE")),
  list(category = "Small_Molecule", reason = "ATC L01EL (Bruton's tyrosine kinase (BTK) inhibitors): targeted small-molecule inhibitor",
       drugs = c("ACALABRUTINIB", "IBRUTINIB")),
  list(category = "Small_Molecule", reason = "ATC L01EX (Other protein kinase inhibitors): targeted small-molecule inhibitor",
       drugs = c("CABOZANTINIB SMALATE", "ENTRECTINIB", "SORAFENIB TOSYLATE", "SUNITINIB MALATE")),
  list(category = "Small_Molecule", reason = "JUDGEMENT CALL: ATC L01XG (Proteasome inhibitors); targeted oral small molecule, though not a kinase inhibitor",
       drugs = c("BORTEZOMIB", "CARFILZOMIB")),
  list(category = "Small_Molecule", reason = "JUDGEMENT CALL: ATC L01XK (PARP inhibitors); targeted oral small molecule, though not a kinase inhibitor",
       drugs = c("OLAPARIB")),
  list(category = "Small_Molecule", reason = "JUDGEMENT CALL: no ATC code (not in RxNorm); targeted small-molecule inhibitor (ensartinib ALK, icotinib + rociletinib EGFR, navitoclax BCL-2)",
       drugs = c("ENSARTINIB", "ICOTINIB", "NAVITOCLAX", "ROCILETINIB")),
  list(category = "Targeted", reason = "ATC L01FA (CD20 (Clusters of Differentiation 20) inhibitors): targeted monoclonal antibody",
       drugs = c("RITUXIMAB", "RITUXIMAB AND HYALURONIDASE HUMAN")),
  list(category = "Targeted", reason = "ATC L01FC (CD38 (Clusters of Differentiation 38) inhibitors): targeted monoclonal antibody",
       drugs = c("DARATUMUMAB")),
  list(category = "Targeted", reason = "ATC L01FD (HER2 (Human Epidermal Growth Factor Receptor 2) inhibitors): targeted monoclonal antibody",
       drugs = c("TRASTUZUMAB", "TRASTUZUMAB EMTANSINE")),
  list(category = "Targeted", reason = "ATC L01FE (EGFR (Epidermal Growth Factor Receptor) inhibitors): targeted monoclonal antibody",
       drugs = c("CETUXIMAB", "NECITUMUMAB")),
  list(category = "Targeted", reason = "ATC L01FG (VEGF/VEGFR (Vascular Endothelial Growth Factor / -Receptor) inhibitors): targeted monoclonal antibody",
       drugs = c("BEVACIZUMAB", "RAMUCIRUMAB")),
  list(category = "IO", reason = "ATC L01FF (PD-1/PD-L1 (Programmed cell death protein 1/death ligand 1) inhibitors): immunotherapy",
       drugs = c("ATEZOLIZUMAB", "AVELUMAB", "DURVALUMAB", "NIVOLUMAB", "PEMBROLIZUMAB")),
  list(category = "IO", reason = "JUDGEMENT CALL: ATC L01FX (Other monoclonal antibodies), but these are CTLA-4 checkpoint inhibitors",
       drugs = c("IPILIMUMAB", "TREMELIMUMAB")),
  list(category = "IO", reason = "JUDGEMENT CALL: ATC L01XL (Antineoplastic cell and gene therapy); oncolytic virus that works through the immune system",
       drugs = c("TALIMOGENE LAHERPAREPVEC")),
  list(category = "IO", reason = "ATC L03AC (Interleukins): immunotherapy",
       drugs = c("ALDESLEUKIN")),
  list(category = "IO", reason = "JUDGEMENT CALL: no ATC link in RxClass; WHO ATC lists BCG as L03AX03 (immunostimulant)",
       drugs = c("BCG SOLUTION", "BCG VACCINE")),
  list(category = "Hormone", reason = "ATC H01CB (Somatostatin and analogues): endocrine therapy",
       drugs = c("LANREOTIDE ACETATE")),
  list(category = "Hormone", reason = "ATC L02AB (Progestogens): endocrine therapy",
       drugs = c("MEGESTROL ACETATE")),
  list(category = "Hormone", reason = "ATC L02AE (Gonadotropin releasing hormone analogues): endocrine therapy",
       drugs = c("LEUPROLIDE ACETATE")),
  list(category = "Hormone", reason = "ATC L02BA (Anti-estrogens): endocrine therapy",
       drugs = c("TAMOXIFEN CITRATE")),
  list(category = "Hormone", reason = "ATC L02BB (Anti-androgens): endocrine therapy",
       drugs = c("BICALUTAMIDE")),
  list(category = "Hormone", reason = "ATC L02BG (Aromatase inhibitors): endocrine therapy",
       drugs = c("ANASTROZOLE", "EXEMESTANE", "LETROZOLE")),
  list(category = "Hormone", reason = "ATC L02BX (Other hormone antagonists and related agents): endocrine therapy",
       drugs = c("ABIRATERONE ACETATE", "DEGARELIX")),
  list(category = "Hormone", reason = "No ATC code (placeholder name); the name says it is a hormone therapy",
       drugs = c("OTHER HORMONE")),
  list(category = "Ancillary", reason = "ATC V03AF (Detoxifying agents for antineoplastic treatment): supportive care, not anticancer",
       drugs = c("LEUCOVORIN CALCIUM")),
  list(category = "Others", reason = "JUDGEMENT CALL: ATC L01XF (Retinoids for cancer treatment); differentiation agent that fits no other category",
       drugs = c("TRETINOIN")),
  list(category = "Others", reason = "JUDGEMENT CALL: ATC L04AX (Other immunosuppressants); immunomodulator that fits no other category",
       drugs = c("LENALIDOMIDE", "POMALIDOMIDE")),
  list(category = "Others", reason = "No ATC code: placeholder name, the actual drug is unknown",
       drugs = c("INVESTIGATIONAL DRUG", "OTHER ANTINEOPLASTIC", "OTHER NOS"))
)

# drug name -> category, and drug name -> reason
all_drugs <- unlist(lapply(groups, function(g) g$drugs))
category_of <- setNames(unlist(lapply(groups, function(g) rep(g$category, length(g$drugs)))), all_drugs)
reason_of <- setNames(unlist(lapply(groups, function(g) rep(g$reason, length(g$drugs)))), all_drugs)
if (any(duplicated(all_drugs))) stop("Drug in more than one group: ", paste(all_drugs[duplicated(all_drugs)], collapse = ", "))

data <- read.csv(file.path(ROOT, "Data", "regimen_cancer_level_dataset.csv"),
                 stringsAsFactors = FALSE, na.strings = c("", "NA"))

# one row per drug: each regimen row has up to 5 drugs (drugs_drug_1 .. drugs_drug_5)
drugs <- do.call(rbind, lapply(1:5, function(n) {
  d <- data[[paste0("drugs_drug_", n)]]
  has_drug <- !is.na(d)
  data.frame(record_id = data$record_id[has_drug],
             drug = toupper(trimws(sub("\\(.*$", "", d[has_drug]))))   # text before "("
}))

drugs$category <- unname(category_of[drugs$drug])
drugs$reason <- unname(reason_of[drugs$drug])

unmatched <- unique(drugs$drug[is.na(drugs$category)])
if (length(unmatched) > 0) stop("No category for: ", paste(unmatched, collapse = ", "))

write.csv(drugs, file.path(ROOT, "Timelines", "R_scripts", "patient_drug_categories.csv"), row.names = FALSE)

cat(nrow(drugs), "drugs from", length(unique(drugs$record_id)), "patients, all matched\n")
print(table(drugs$category))
