"""Step 09 -- Build the audited TCGA-COAD cohort from raw clinical + expression files.

This step was previously missing from the repository -- it existed only as an
interactive analysis during development and was never saved as a runnable script.
Reconstructed here so the full pipeline is actually reproducible end to end.

Expects, per docs/data_provenance.md:
    data/raw/TCGA-COAD.clinical.tsv.gz
    data/raw/TCGA-COAD.star_counts.tsv.gz

Produces:
    tcga/tumor_final_audit.csv   (454 patients, one Primary Tumor sample each, valid age)
"""
import os
import pandas as pd

os.makedirs("tcga", exist_ok=True)

clinical = pd.read_csv("data/raw/TCGA-COAD.clinical.tsv.gz", sep="\t", low_memory=False)
counts_header = pd.read_csv("data/raw/TCGA-COAD.star_counts.tsv.gz", sep="\t", nrows=0)
expr_samples = set(counts_header.columns[1:])  # exclude the Ensembl_ID column

sub = clinical[clinical["sample"].isin(expr_samples)].copy()
sub["patient_id"] = sub["sample"].str[:12]
tumor = sub[sub["sample_type.samples"] == "Primary Tumor"].copy()
print(f"Samples with both clinical + expression data: {len(sub)}")
print(f"Primary Tumor samples: {len(tumor)}, unique patients: {tumor.patient_id.nunique()}")

# 13 patients have 2 Primary Tumor samples -- keep the most recently collected one
# (documented, arbitrary-but-consistent tie-break; see docs/decision_log.md)
tumor = tumor.sort_values("days_to_collection.samples", ascending=False) \
             .drop_duplicates("patient_id", keep="first")

age_years = pd.to_numeric(tumor["age_at_earliest_diagnosis_in_years.diagnoses.xena_derived"],
                           errors="coerce")
tumor["age_years"] = age_years
tumor["age_group"] = pd.NA
tumor.loc[age_years < 50, "age_group"] = "EOCRC"
tumor.loc[age_years >= 50, "age_group"] = "LOCRC"

tumor_valid = tumor.dropna(subset=["age_group"]).copy()
print(f"\nFinal audited TCGA cohort: {len(tumor_valid)} patients "
      f"({(tumor_valid.age_group=='EOCRC').sum()} EOCRC, "
      f"{(tumor_valid.age_group=='LOCRC').sum()} LOCRC)")

tumor_valid.to_csv("tcga/tumor_final_audit.csv", index=False)
print("Saved tcga/tumor_final_audit.csv")
