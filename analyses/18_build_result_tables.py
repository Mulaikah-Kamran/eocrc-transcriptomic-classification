"""Step 18 -- Build the four core summary tables referenced throughout REPORT.md.

Table 4 (full Hallmark enrichment) is produced directly by 16_pathway_enrichment.py.
This script builds Tables 1-3.

Reads:
    data/primary_metadata.csv, tcga/tumor_final_audit.csv
    data/repA_shap_importance.csv, data/repA_gene_stability.csv
    data/tcga_prediction_scores_representation_a.csv (etc., for AUC values -- see note below)
Produces:
    results/tables/table1_cohort_characteristics.csv
    results/tables/table2_representation_comparison.csv
    results/tables/table3_top_stable_genes.csv
"""
import os
import pandas as pd

os.makedirs("results/tables", exist_ok=True)

# --- Table 1: cohort characteristics ---
meta = pd.read_csv("data/primary_metadata.csv")
tcga = pd.read_csv("tcga/tumor_final_audit.csv")

table1 = pd.DataFrame({
    "Cohort": ["GSE213092 (discovery)", "TCGA-COAD (external validation)"],
    "n": [len(meta), len(tcga)],
    "EOCRC": [(meta.cohort == "EOCRC").sum(), (tcga.age_group == "EOCRC").sum()],
    "LOCRC": [(meta.cohort == "LOCRC").sum(), (tcga.age_group == "LOCRC").sum()],
    "Age range (EOCRC)": ["25-49", "<50"],
    "Age range (LOCRC)": ["71-80", "31-90 (>=50)"],
    "Platform": ["Illumina HiSeq 2500 (single)", "GDC-harmonized (multi-site, STAR pipeline)"],
    "Notes": ["Single BioProject, tumor-only",
              "Severe class imbalance; site-vs-age p=0.0045 (partial, non-fatal confound)"],
})
table1.to_csv("results/tables/table1_cohort_characteristics.csv", index=False)

# --- Table 2: representation comparison ---
# NOTE: these AUC/p-value figures are the frozen, reported results from this project's
# full analysis run (see REPORT.md Section 4.1). Re-running the pipeline end-to-end
# reproduces them; they are not re-computed live in this table-building step, since
# doing so would require re-reading results from three separate prediction files
# with three different permutation procedures -- kept explicit here instead.
table2 = pd.DataFrame({
    "Representation": ["Gene-level (2,000 genes)", "Hallmark pathway (50 sets)", "Marx 8-gene panel"],
    "Feature construction": [
        "Robust-expression filter + top-2000 variance genes, unsupervised",
        "Mean z-score of Hallmark member genes per pathway",
        "8 externally nominated genes, no filtering/selection",
    ],
    "Model": ["Elastic Net logistic regression (nested CV)",
              "Elastic Net logistic regression (nested CV)",
              "Ridge logistic regression (nested CV)"],
    "Discovery CV AUC": [0.666, 0.624, 0.591],
    "External TCGA AUC": [0.626, 0.548, 0.540],
    "External permutation p": [0.001, 0.127, 0.174],
})
table2.to_csv("results/tables/table2_representation_comparison.csv", index=False)

# --- Table 3: top 30 stable genes by SHAP importance ---
shap_imp = pd.read_csv("data/repA_shap_importance.csv", index_col=0).iloc[:, 0].sort_values(ascending=False)
stability = pd.read_csv("data/repA_gene_stability.csv", index_col=0)
table3 = pd.DataFrame({"gene": shap_imp.head(30).index, "shap_importance": shap_imp.head(30).values})
table3["nonzero_rate_bootstrap"] = table3["gene"].map(stability["nonzero_rate"])
table3["selection_rate_bootstrap"] = table3["gene"].map(stability["selection_rate"])
table3.to_csv("results/tables/table3_top_stable_genes.csv", index=False)

print("Saved table1, table2, table3 to results/tables/")
