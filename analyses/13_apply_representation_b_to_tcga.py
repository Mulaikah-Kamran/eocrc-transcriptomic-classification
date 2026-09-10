"""Step 13 -- Apply the frozen Representation B model to TCGA-COAD.

Reuses tcga/tcga_log2cpm_symbols.csv from 11_build_tcga_log2cpm.py.

Reads:
    data/frozen_representation_b_model.pkl
    tcga/tcga_log2cpm_symbols.csv
    tcga/tumor_final_audit.csv
Produces:
    data/tcga_prediction_scores_representation_b.csv
"""
import pandas as pd
import numpy as np
import pickle
from sklearn.metrics import (roc_auc_score, average_precision_score, balanced_accuracy_score,
                              confusion_matrix, recall_score)

with open("data/frozen_representation_b_model.pkl", "rb") as f:
    frozen_b = pickle.load(f)

log2cpm_tcga = pd.read_csv("tcga/tcga_log2cpm_symbols.csv", index_col=0)
tumor_valid = pd.read_csv("tcga/tumor_final_audit.csv")

gene_mean = frozen_b["gene_mean"]
gene_std = frozen_b["gene_std"]

common_genes = [g for g in gene_mean.index if g in log2cpm_tcga.index]
z_tcga = log2cpm_tcga.loc[common_genes].sub(gene_mean[common_genes], axis=0).div(gene_std[common_genes], axis=0)

# pathway score = mean z-score of AVAILABLE member genes only -- missing genes are
# excluded from the average rather than treated as zero (see docs/decision_log.md
# for why this differs from Representation A's missing-gene convention)
pathway_scores_tcga = {}
for name, genes in frozen_b["genesets"].items():
    avail = [g for g in genes if g in z_tcga.index]
    pathway_scores_tcga[name] = z_tcga.loc[avail].mean(axis=0)
X_tcga = pd.DataFrame(pathway_scores_tcga)[frozen_b["pathway_names"]]

X_tcga_s = (X_tcga.values - frozen_b["pathway_scaler_mean"]) / frozen_b["pathway_scaler_scale"]

y_true = tumor_valid.set_index("sample").loc[X_tcga.index, "age_group"].map({"EOCRC": 1, "LOCRC": 0}).values
proba = frozen_b["model"].predict_proba(X_tcga_s)[:, 1]
pred = (proba >= 0.5).astype(int)

auc = roc_auc_score(y_true, proba)
pr_auc = average_precision_score(y_true, proba)
bal_acc = balanced_accuracy_score(y_true, pred)
sens = recall_score(y_true, pred, pos_label=1)
spec = recall_score(y_true, pred, pos_label=0)
cm = confusion_matrix(y_true, pred)

print(f"=== Representation B: TCGA-COAD external validation ===")
print(f"Patients evaluated: {len(y_true)}")
print(f"ROC-AUC (primary): {auc:.3f} | PR-AUC: {pr_auc:.3f} | Balanced accuracy: {bal_acc:.3f}")
print(f"Sensitivity: {sens:.3f} | Specificity: {spec:.3f}")
print(f"Confusion matrix:\n{cm}")

pd.DataFrame({"proba": proba, "true": np.where(y_true == 1, "EOCRC", "LOCRC")}) \
    .to_csv("data/tcga_prediction_scores_representation_b.csv", index=False)
