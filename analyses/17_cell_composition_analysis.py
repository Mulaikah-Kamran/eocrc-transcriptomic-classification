"""Step 17 -- Pre-specified canonical cell-type marker composition analysis.

Marker panel and all three tests were decided and documented BEFORE examining
any of these results -- see docs/decision_log.md.

Reads:
    data/primary_log2cpm.csv, data/primary_fold_assignment.csv
    data/repA_high_confidence_genes.csv, data/frozen_representation_a_model.pkl
Produces:
    results/tables/cell_composition_results.csv
    data/cell_composition_scores.csv
"""
import os
import pandas as pd
import numpy as np
import pickle
from scipy.stats import mannwhitneyu, spearmanr

os.makedirs("results/tables", exist_ok=True)

MARKERS = {
    "Immune/leukocyte": ["PTPRC"],
    "Fibroblast/stromal": ["FAP", "PDGFRB"],
    "Endothelial/vascular": ["PECAM1", "VWF"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18"],
    "Smooth muscle/pericyte": ["ACTA2", "MYH11", "TAGLN"],  # overlaps myofibroblast markers -- see docs/decision_log.md
}

log2cpm = pd.read_csv("data/primary_log2cpm.csv", index_col=0)
folds = pd.read_csv("data/primary_fold_assignment.csv", index_col=0)
high_conf = pd.read_csv("data/repA_high_confidence_genes.csv", index_col=0)

with open("data/frozen_representation_a_model.pkl", "rb") as f:
    frozen_a = pickle.load(f)
all_patients = folds.index.tolist()
X_all = log2cpm.loc[frozen_a["final_gene_list"], all_patients].T.values
X_all_s = (X_all - frozen_a["scaler_mean"]) / frozen_a["scaler_scale"]
classifier_proba = pd.Series(frozen_a["model"].predict_proba(X_all_s)[:, 1], index=all_patients)

gene_mean = log2cpm[all_patients].mean(axis=1)
gene_std = log2cpm[all_patients].std(axis=1).replace(0, np.nan)
z_all = log2cpm[all_patients].sub(gene_mean, axis=0).div(gene_std, axis=0)

cohort = folds["cohort"]
composite_scores, results = {}, []
for compartment, genes in MARKERS.items():
    present_in_722 = [g for g in genes if g in high_conf.index]
    avail = [g for g in genes if g in z_all.index]
    score = z_all.loc[avail, all_patients].mean(axis=0)
    composite_scores[compartment] = score
    g1, g2 = score[cohort == "EOCRC"], score[cohort == "LOCRC"]
    _, p_group = mannwhitneyu(g1, g2)
    rho, p_corr = spearmanr(classifier_proba.loc[all_patients], score.loc[all_patients])
    results.append(dict(compartment=compartment, markers_present_in_high_confidence_set=present_in_722,
                         eocrc_mean=g1.mean(), locrc_mean=g2.mean(), group_diff_p=p_group,
                         spearman_rho_vs_classifier=rho, spearman_p=p_corr))
    print(f"{compartment}: present in high-confidence set = {present_in_722}, "
          f"group-diff p={p_group:.3f}, corr-with-model p={p_corr:.3f}")

pd.DataFrame(results).to_csv("results/tables/cell_composition_results.csv", index=False)
pd.DataFrame(composite_scores).to_csv("data/cell_composition_scores.csv")
print("\nSaved results/tables/cell_composition_results.csv and data/cell_composition_scores.csv")
