"""Step 12 -- Apply the frozen Representation A model to TCGA-COAD.

Reuses tcga/tcga_log2cpm_symbols.csv from 11_build_tcga_log2cpm.py rather than
re-deriving it (fixes duplication from an earlier version of this pipeline).

No retraining, no retuning, no feature selection using TCGA data at any point.

Reads:
    data/frozen_representation_a_model.pkl
    tcga/tcga_log2cpm_symbols.csv
    tcga/tumor_final_audit.csv
Produces:
    data/tcga_prediction_scores_representation_a.csv
"""
import pandas as pd
import numpy as np
import pickle
from sklearn.metrics import (roc_auc_score, average_precision_score, balanced_accuracy_score,
                              confusion_matrix, recall_score)

with open("data/frozen_representation_a_model.pkl", "rb") as f:
    frozen = pickle.load(f)
final_genes = frozen["final_gene_list"]
scaler_mean = pd.Series(frozen["scaler_mean"], index=final_genes)
scaler_scale = pd.Series(frozen["scaler_scale"], index=final_genes)

log2cpm_tcga = pd.read_csv("tcga/tcga_log2cpm_symbols.csv", index_col=0)
tumor_valid = pd.read_csv("tcga/tumor_final_audit.csv")

X_tcga = pd.DataFrame(index=tumor_valid["sample"], columns=final_genes, dtype=float)
available_genes = [g for g in final_genes if g in log2cpm_tcga.index]
missing_genes = [g for g in final_genes if g not in log2cpm_tcga.index]
print(f"Genes available in TCGA: {len(available_genes)} / {len(final_genes)}")

for g in available_genes:
    X_tcga[g] = log2cpm_tcga.loc[g, tumor_valid["sample"]].values

X_tcga_z = (X_tcga - scaler_mean) / scaler_scale
# documented convention: missing genes contribute nothing (standardized value = 0)
X_tcga_z[missing_genes] = 0.0

y_true = tumor_valid.set_index("sample").loc[X_tcga_z.index, "age_group"].map({"EOCRC": 1, "LOCRC": 0}).values
proba = frozen["model"].predict_proba(X_tcga_z[final_genes].values)[:, 1]
pred = (proba >= 0.5).astype(int)

auc = roc_auc_score(y_true, proba)
pr_auc = average_precision_score(y_true, proba)
bal_acc = balanced_accuracy_score(y_true, pred)
sens = recall_score(y_true, pred, pos_label=1)
spec = recall_score(y_true, pred, pos_label=0)
cm = confusion_matrix(y_true, pred)
majority_baseline = max((y_true == 1).mean(), (y_true == 0).mean())

print(f"\n=== Representation A: TCGA-COAD external validation ===")
print(f"Patients evaluated: {len(y_true)} (EOCRC={int((y_true==1).sum())}, LOCRC={int((y_true==0).sum())})")
print(f"ROC-AUC (primary): {auc:.3f} | PR-AUC: {pr_auc:.3f} | Balanced accuracy: {bal_acc:.3f}")
print(f"Sensitivity: {sens:.3f} | Specificity: {spec:.3f}")
print(f"Confusion matrix:\n{cm}")
print(f"Majority-class baseline: {majority_baseline:.3f}")

pd.DataFrame({"proba": proba, "true": np.where(y_true == 1, "EOCRC", "LOCRC")}) \
    .to_csv("data/tcga_prediction_scores_representation_a.csv", index=False)
