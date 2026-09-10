"""Step 08 -- Freeze the final Representation B (Hallmark pathway) model on all 98 patients.

Reads:
    data/primary_log2cpm.csv, data/primary_fold_assignment.csv
    data/raw/hallmark_gene_sets.gmt
Produces: data/frozen_representation_b_model.pkl
"""
import sys
import pandas as pd
import numpy as np
import pickle
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, ".")
from src.gmt_utils import parse_gmt

log2cpm = pd.read_csv("data/primary_log2cpm.csv", index_col=0)
folds = pd.read_csv("data/primary_fold_assignment.csv", index_col=0)
genesets = parse_gmt("data/raw/hallmark_gene_sets.gmt", restrict_to=log2cpm.index)

y_map = {"EOCRC": 1, "LOCRC": 0}
all_patients = folds.index.tolist()
y_all = folds["cohort"].map(y_map).values

gene_mean = log2cpm[all_patients].mean(axis=1)
gene_std = log2cpm[all_patients].std(axis=1).replace(0, np.nan)
z_all = log2cpm[all_patients].sub(gene_mean, axis=0).div(gene_std, axis=0)

pathway_scores = pd.DataFrame({name: z_all.loc[list(genes)].mean(axis=0) for name, genes in genesets.items()})
X_all = pathway_scores.loc[all_patients].values
print(f"Pathway score matrix: {X_all.shape[0]} patients x {X_all.shape[1]} pathways")

scaler = StandardScaler().fit(X_all)
X_all_s = scaler.transform(X_all)

inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
final_model = LogisticRegressionCV(
    penalty="elasticnet", solver="saga",
    l1_ratios=[0.3, 0.7], Cs=5,
    cv=inner_cv, max_iter=1000, tol=1e-3, scoring="roc_auc",
    random_state=42, n_jobs=1
)
final_model.fit(X_all_s, y_all)

print(f"Chosen hyperparameters: C={final_model.C_[0]:.4f}, l1_ratio={final_model.l1_ratio_[0]}")
print(f"Pathways with non-zero coefficient: {np.sum(final_model.coef_[0] != 0)} / 50")

with open("data/frozen_representation_b_model.pkl", "wb") as f:
    pickle.dump({
        "genesets": genesets,
        "gene_mean": gene_mean,
        "gene_std": gene_std,
        "pathway_scaler_mean": scaler.mean_,
        "pathway_scaler_scale": scaler.scale_,
        "pathway_names": list(pathway_scores.columns),
        "model": final_model,
    }, f)
print("Saved data/frozen_representation_b_model.pkl")
