# Runs the MEC-TX R package (mectx, CRAN) on the GENIE NSCLC regimen data.
# Run from anywhere:   Rscript mectx_run.R   (or Source in VS Code / RStudio)
# Outputs: mectx_clusters.csv in Timelines/R_scripts, plots in Timelines/R_scripts/plots.
#
# mectx expects ORIEN-style columns, so this script converts the GENIE CSV first:
#   sample            <- record_id
#   AgeAtMedStart/Stop <- drug start/end, in YEARS SINCE DIAGNOSIS (mectx only uses differences, so
#                        "years since diagnosis" works the same as age)
#   AgeAtLastContact  <- death / last contact = regimen drug-1 start + tt_os_d1_days
#   treatment_group   <- drug class (Chemo, IO, Targeted, Hormone, Ancillary, Others)
#   diagsurvtime      <- overall survival from diagnosis, in years;  status <- os_d_status
suppressPackageStartupMessages({
  library(mectx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ALEX - This is just for me when im remote vs at home
pathlaptop <- "C:/Users/default.DESKTOP-UOBTTU9/Bioinformatics Lesgo/rxAlign"
#pathhome <- "C:/Users/Zenith/PhDCoding/rxAlign"

ROOT <- pathlaptop
#ROOT <- pathhome

#setting up the paths based on root
DATA_PATH <- file.path(ROOT, "data", "regimen_cancer_level_dataset.csv")
MED_CLASS_PATH <- file.path(ROOT, "Code", "data", "med_class.tsv")
OUT_DIR <- file.path(ROOT, "Timelines", "R_scripts")       # CSV output
PLOT_DIR <- file.path(OUT_DIR, "plots")                     # plot output

dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

FINAL_K <- 5            # k to plot / run survival on (all k from 3 to 20 are saved)
INDEX_CANCER_ONLY <- TRUE    # skip regimens for a patient's other cancers (redcap_ca_index == "No")
DAYS_PER_YEAR <- 365.25 #needed for age calculations, is the .25 fine?

# SETUP TO MECTX FORMAT:
# 1. drug -> treatment class 
# start from byu group med_class.tsv (Targeted / Chemo / ICB / Other)
med_class <- read.delim(MED_CLASS_PATH, stringsAsFactors = FALSE)
drug_class <- setNames(
  c(Targeted = "Targeted", Chemo = "Chemo", ICB = "IO", Other = "Others")[med_class$Class],
  toupper(trimws(med_class$Medication))
)
# drugs missing from med_class.tsv, plus hormone therapies (mectx has a separate Hormone class)
overrides <- list(
  Targeted = c("AFATINIB DIMALEATE", "CABOZANTINIB SMALATE", "DARATUMUMAB", "ENSARTINIB", "ENTRECTINIB",
               "ERLOTINIB HYDROCHLORIDE", "ICOTINIB", "NAVITOCLAX", "NECITUMUMAB", "RITUXIMAB AND HYALURONIDASE HUMAN",
               "ROCILETINIB", "RUXOLITINIB PHOSPHATE", "TRASTUZUMAB EMTANSINE"),
  Chemo = c("ARSENIC TRIOXIDE", "BLEOMYCIN", "BUSULFAN", "CHLORAMBUCIL", "ETOPOSIDE PHOSPHATE", "LURBINECTEDIN",
            "MECHLORETHAMINE HYDROCHLORIDE", "PACLITAXEL LOADED POLYMERIC MICELLE", "PACLITAXEL POLIGLUMEX",
            "PACLITAXEL TREVATIDE", "PEGYLATED LIPOSOMAL DOXORUBICIN HYDROCHLORIDE", "VINCRISTINE SULFATE LIPOSOME"),
  IO = c("TREMELIMUMAB", "BCG SOLUTION", "BCG VACCINE"),
  Hormone = c("ANASTROZOLE", "LETROZOLE", "TAMOXIFEN CITRATE", "EXEMESTANE", "LEUPROLIDE ACETATE", "BICALUTAMIDE",
              "DEGARELIX", "MEGESTROL ACETATE", "ABIRATERONE ACETATE", "LANREOTIDE ACETATE", "OTHER HORMONE"),
  Ancillary = c("LEUCOVORIN CALCIUM"),
  Others = c("INVESTIGATIONAL DRUG", "OTHER ANTINEOPLASTIC", "OTHER NOS")
)
for (cls in names(overrides)) drug_class[overrides[[cls]]] <- cls

# ---------- read + reshape GENIE data ----------
raw <- read.csv(DATA_PATH, stringsAsFactors = FALSE, na.strings = c("", "NA"))
if (INDEX_CANCER_ONLY) raw <- raw %>% filter(redcap_ca_index == "Yes")

# one row per drug (drugs_drug_1..5 -> long)
drugs <- bind_rows(lapply(1:5, function(n) {
  raw %>% transmute(
    sample = record_id,
    drug = .data[[paste0("drugs_drug_", n)]],
    start = .data[[paste0("dx_drug_start_int_", n)]],
    # end date; if blank, fall back to last administration (if still blank, mectx uses last contact)
    end = coalesce(.data[[paste0("dx_drug_end_int_", n)]], .data[[paste0("dx_drug_end_or_lastadm_int_", n)]]),
    last_contact = dx_drug_start_int_1 + tt_os_d1_days
  )
})) %>% filter(!is.na(drug), !is.na(start))

base_name <- toupper(trimws(sub("\\(.*$", "", drugs$drug)))   # drop the "(synonyms...)"
unmapped <- setdiff(unique(base_name), names(drug_class))
if (length(unmapped) > 0) message("No class for (using Others): ", paste(unmapped, collapse = ", "))

med_data <- drugs %>%
  mutate(Medication = base_name,
         treatment_group = unname(ifelse(base_name %in% names(drug_class), drug_class[base_name], "Others"))) %>%
  group_by(sample) %>%
  mutate(AgeAtLastContact = max(last_contact, na.rm = TRUE) / DAYS_PER_YEAR) %>%
  ungroup() %>%
  transmute(sample, Medication, treatment_group,
            AgeAtMedStart = start / DAYS_PER_YEAR,
            AgeAtMedStop = end / DAYS_PER_YEAR,
            Age.At.Specimen.Collection = NA_real_,   # not used here: don't clip treatment to specimen date
            AgeAtLastContact)

metadata <- raw %>%
  group_by(sample = record_id) %>%
  summarise(diagsurvtime = max(dx_drug_start_int_1 + tt_os_d1_days, na.rm = TRUE) / DAYS_PER_YEAR,
            status = first(os_d_status), .groups = "drop")

message(sprintf("%d patients, %d drug records", n_distinct(med_data$sample), nrow(med_data)))

# ---------- MEC-TX pipeline ----------
norm <- tx_normalize(med_data, metadata)
res <- tx_cluster_surv(metadata, norm)
write.csv(res$Cluster_surv, file.path(OUT_DIR, "mectx_clusters.csv"), row.names = FALSE)

# silhouette for every k (mectx returns all k but doesn't pick one)
d <- dist(res$pca_matrix)
k_cols <- grep("^Cluster_k", names(res$cluster_results), value = TRUE)
sil <- sapply(k_cols, function(kc) mean(cluster::silhouette(res$cluster_results[[kc]], d)[, 3]))
print(data.frame(k = sub("Cluster_k", "", k_cols), silhouette = round(sil, 3)), row.names = FALSE)

kc <- paste0("Cluster_k", FINAL_K)

# UMAP colored by cluster
umap_plot <- res$umap_df %>%
  left_join(res$cluster_results %>% select(sample, cluster = all_of(kc)), by = "sample") %>%
  ggplot(aes(UMAP1, UMAP2, color = factor(cluster))) +
  geom_point(size = 0.8, alpha = 0.6) +
  labs(color = "cluster", title = paste("MEC-TX UMAP,", kc)) +
  theme_bw()
ggsave(file.path(PLOT_DIR, "mectx_umap.png"), umap_plot, width = 7, height = 5, dpi = 120)

# swimlane timelines per cluster
segs <- tx_intervals(norm)
timeline_plot <- plot_timeline_for_k(kc, res$cluster_results, segs, horizon_years = 5, base_size = 10, title_size = 14)
ggsave(file.path(PLOT_DIR, "mectx_timelines.png"), timeline_plot, width = 14, height = 8, dpi = 120)

# survival by cluster (Kaplan-Meier + Cox); no stage/sex/age columns in this file, so Cox is cluster-only
surv <- tx_compare_groups(res$Cluster_surv, group_var = kc)
ggsave(file.path(PLOT_DIR, "mectx_survival.png"), surv$combined, width = 14, height = 6, dpi = 120)
print(surv$group_summary)
