"""Step 05 -- Initial (40-shuffle) permutation test on the discovery-cohort result.

This is the first, smaller permutation test -- see docs/decision_log.md and
docs/limitations.md for why this was later extended to 1,200 shuffles
(06_permutation_test_discovery_1200.py) rather than trusted on its own.

Uses a FIXED-hyperparameter version of the Representation A pipeline (no inner
tuning loop) purely for speed -- this is a secondary significance check, not
the frozen model itself.

Reads: data/primary_counts_raw.csv, data/primary_log2cpm.csv, data/primary_fold_assignment.csv
Produces: data/permutation_discovery_40shuffle.csv
"""
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings("ignore")

log2cpm = pd.read_csv("data/primary_log2cpm.csv", index_col=0)
folds = pd.read_csv("data/primary_fold_assignment.csv", index_col=0)
counts = pd.read_csv("data/primary_counts_raw.csv", index_col=0)
libsize = counts.sum(axis=0)
cpm = counts.div(libsize, axis=1) * 1e6
y_map = {"EOCRC": 1, "LOCRC": 0}
N_TOP_GENES = 2000
true_labels = folds.cohort.map(y_map)


def run_cv_fixed_hparams(label_series, C=1.0, l1_ratio=0.5):
    all_true, all_proba = [], []
    for fold_i in sorted(folds.fold.unique()):
        train_p = folds[folds.fold != fold_i].index.tolist()
        test_p = folds[folds.fold == fold_i].index.tolist()
        train_cohort = label_series.loc[train_p]
        eocrc_train = train_cohort[train_cohort == 1].index
        locrc_train = train_cohort[train_cohort == 0].index
        detected_e = (cpm.loc[:, eocrc_train] > 1.0).mean(axis=1) >= 0.7
        detected_l = (cpm.loc[:, locrc_train] > 1.0).mean(axis=1) >= 0.7
        robust_genes = cpm.index[detected_e & detected_l]
        top_genes = log2cpm.loc[robust_genes, train_p].var(axis=1).nlargest(N_TOP_GENES).index

        X_train = log2cpm.loc[top_genes, train_p].T.values
        X_test = log2cpm.loc[top_genes, test_p].T.values
        y_train = train_cohort.values
        y_test = label_series.loc[test_p].values

        scaler = StandardScaler().fit(X_train)
        clf = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=l1_ratio,
                                  C=C, max_iter=500, tol=1e-2, random_state=0)
        clf.fit(scaler.transform(X_train), y_train)
        proba = clf.predict_proba(scaler.transform(X_test))[:, 1]
        all_true.extend(y_test)
        all_proba.extend(proba)
    return roc_auc_score(all_true, all_proba)


real_auc = run_cv_fixed_hparams(true_labels)
print(f"Real (unpermuted) fixed-hyperparameter AUC: {real_auc:.4f}")

N_PERM = 40
rng = np.random.default_rng(42)
perm_aucs = []
for i in range(N_PERM):
    shuffled = pd.Series(rng.permutation(true_labels.values), index=true_labels.index)
    perm_aucs.append(run_cv_fixed_hparams(shuffled))
    if (i + 1) % 10 == 0:
        print(f"  {i+1}/{N_PERM} permutations done")

perm_aucs = np.array(perm_aucs)
p_val = (np.sum(perm_aucs >= real_auc) + 1) / (len(perm_aucs) + 1)
print(f"\nPermutation null AUCs: mean={perm_aucs.mean():.3f}, max={perm_aucs.max():.3f}")
print(f"Empirical p-value (n={N_PERM}): {p_val:.4f}")

pd.Series(perm_aucs, name="perm_auc").to_csv("data/permutation_discovery_40shuffle.csv", index=False)
