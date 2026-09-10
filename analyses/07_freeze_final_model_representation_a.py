"""Step 07 -- Freeze the final Representation A model on all 98 discovery patients.

This is the ONE model that gets applied to TCGA later -- no refitting, no retuning
after this point. See docs/decision_log.md for why this is a separate step from
the cross-validation in 03 (that estimates performance; this produces the actual
model used going forward).

Reads: data/primary_counts_raw.csv, data/primary_log2cpm.csv, data/primary_fold_assignment.csv
Produces: data/frozen_representation_a_model.pkl
"""
import pandas as pd
import numpy as np
import pickle
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
import warnings
warnings.filterwarnings("ignore")

counts = pd.read_csv("data/primary_counts_raw.csv", index_col=0)
log2cpm = pd.read_csv("data/primary_log2cpm.csv", index_col=0)
folds = pd.read_csv("data/primary_fold_assignment.csv", index_col=0)

y_map = {"EOCRC": 1, "LOCRC": 0}
N_TOP_GENES = 2000
MIN_CPM, MIN_FRAC_DETECTED = 1.0, 0.7

all_patients = folds.index.tolist()
y_all = folds["cohort"].map(y_map).values

libsize = counts.sum(axis=0)
cpm = counts.div(libsize, axis=1) * 1e6

eocrc_p = folds[folds.cohort == "EOCRC"].index
locrc_p = folds[folds.cohort == "LOCRC"].index

detected_eocrc = (cpm.loc[:, eocrc_p] > MIN_CPM).mean(axis=1) >= MIN_FRAC_DETECTED
detected_locrc = (cpm.loc[:, locrc_p] > MIN_CPM).mean(axis=1) >= MIN_FRAC_DETECTED
robust_genes = cpm.index[detected_eocrc & detected_locrc]
print(f"Expression-robustness filter: {len(robust_genes)} genes retained")

top_genes = log2cpm.loc[robust_genes, all_patients].var(axis=1).nlargest(N_TOP_GENES).index
print(f"Top-{N_TOP_GENES}-variance selection: {len(top_genes)} genes (final feature set)")

X_all = log2cpm.loc[top_genes, all_patients].T.values
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

print(f"\nFinal model trained on all {len(all_patients)} GSE213092 patients.")
print(f"Chosen hyperparameters: C={final_model.C_[0]:.4f}, l1_ratio={final_model.l1_ratio_[0]}")
n_nonzero = np.sum(final_model.coef_[0] != 0)
print(f"Genes with non-zero coefficient: {n_nonzero} / {len(top_genes)}")

with open("data/frozen_representation_a_model.pkl", "wb") as f:
    pickle.dump({
        "final_gene_list": list(top_genes),
        "scaler_mean": scaler.mean_,
        "scaler_scale": scaler.scale_,
        "model": final_model,
    }, f)
print("Saved data/frozen_representation_a_model.pkl")
