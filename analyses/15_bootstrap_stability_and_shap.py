"""Step 15 -- SHAP importance and bootstrap stability for the frozen Representation A model.

Fixes from the earlier version of this script: removed a dead file-existence
check that opened a file only to immediately discard it (would raise an
unhandled error if the file were ever absent), and this script now produces
data/repA_high_confidence_genes.csv directly, which scripts 16 and 17 depend on
(previously missing).

Reads: data/primary_log2cpm.csv, data/primary_fold_assignment.csv,
       data/primary_counts_raw.csv, data/frozen_representation_a_model.pkl
Produces: data/repA_shap_importance.csv, data/repA_gene_stability.csv,
          data/repA_high_confidence_genes.csv
"""
import pandas as pd
import numpy as np
import pickle
import shap
import time
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")

with open("data/frozen_representation_a_model.pkl", "rb") as f:
    frozen = pickle.load(f)
final_genes = frozen["final_gene_list"]
model = frozen["model"]
scaler_mean = frozen["scaler_mean"]
scaler_scale = frozen["scaler_scale"]

log2cpm = pd.read_csv("data/primary_log2cpm.csv", index_col=0)
folds = pd.read_csv("data/primary_fold_assignment.csv", index_col=0)
counts = pd.read_csv("data/primary_counts_raw.csv", index_col=0)
all_patients = folds.index.tolist()

# --- SHAP (exact, since the frozen model is linear) ---
X_all = log2cpm.loc[final_genes, all_patients].T.values
X_all_s = (X_all - scaler_mean) / scaler_scale
explainer = shap.LinearExplainer(model, X_all_s)
shap_values = explainer(X_all_s)
mean_abs_shap = np.abs(shap_values.values).mean(axis=0)
importance = pd.Series(mean_abs_shap, index=final_genes).sort_values(ascending=False)
importance.to_csv("data/repA_shap_importance.csv")
print(f"SHAP importance computed for {len(importance)} genes")

# --- bootstrap stability (200 resamples, fixed frozen hyperparameters) ---
libsize = counts.sum(axis=0)
cpm = counts.div(libsize, axis=1) * 1e6
eocrc_p = folds[folds.cohort == "EOCRC"].index.tolist()
locrc_p = folds[folds.cohort == "LOCRC"].index.tolist()
original_genes = set(final_genes)

N_BOOT = 200
selection_count = pd.Series(0, index=list(original_genes))
nonzero_count = pd.Series(0, index=list(original_genes))

rng = np.random.default_rng(42)
t0 = time.time()
for b in range(N_BOOT):
    boot_e = list(rng.choice(eocrc_p, size=len(eocrc_p), replace=True))
    boot_l = list(rng.choice(locrc_p, size=len(locrc_p), replace=True))
    idx = boot_e + boot_l
    y = np.array([1] * len(boot_e) + [0] * len(boot_l))

    detected_e = (cpm[boot_e] > 1.0).mean(axis=1) >= 0.7
    detected_l = (cpm[boot_l] > 1.0).mean(axis=1) >= 0.7
    robust = cpm.index[detected_e & detected_l]
    top_genes = log2cpm.loc[robust, idx].var(axis=1).nlargest(2000).index

    selected_here = set(top_genes) & original_genes
    selection_count[list(selected_here)] += 1

    X = log2cpm.loc[top_genes, idx].T.values
    scaler = StandardScaler().fit(X)
    clf = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.3, C=100.0,
                              max_iter=500, tol=1e-2, random_state=0)
    clf.fit(scaler.transform(X), y)
    coefs = pd.Series(clf.coef_[0], index=top_genes)
    nonzero_here = set(coefs[coefs != 0].index) & original_genes
    nonzero_count[list(nonzero_here)] += 1

print(f"{N_BOOT} bootstraps completed in {time.time()-t0:.1f}s")

stability = pd.DataFrame({
    "selection_rate": selection_count / N_BOOT,
    "nonzero_rate": nonzero_count / N_BOOT,
})
stability["shap_importance"] = importance.reindex(stability.index)
stability = stability.sort_values("nonzero_rate", ascending=False)
stability.to_csv("data/repA_gene_stability.csv")

# the "high-confidence" set used throughout the rest of the pipeline (enrichment,
# cell-composition audit): genes with a non-zero coefficient in every single
# bootstrap resample, ranked by SHAP importance
high_confidence = stability[stability.nonzero_rate == 1.0].sort_values("shap_importance", ascending=False)
high_confidence.to_csv("data/repA_high_confidence_genes.csv")
print(f"\nHigh-confidence gene set (nonzero in all {N_BOOT} bootstraps): {len(high_confidence)} genes")
