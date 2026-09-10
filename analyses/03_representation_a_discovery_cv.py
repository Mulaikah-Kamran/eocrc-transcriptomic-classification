"""Step 03 -- Representation A (gene-level), discovery-cohort cross-validation.

Reads: data/primary_counts_raw.csv, data/primary_log2cpm.csv, data/primary_fold_assignment.csv
Produces: data/repA_per_fold_results.csv, data/repA_oof_predictions.csv
"""
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score
import warnings
warnings.filterwarnings("ignore")

counts = pd.read_csv("data/primary_counts_raw.csv", index_col=0)
log2cpm = pd.read_csv("data/primary_log2cpm.csv", index_col=0)
folds = pd.read_csv("data/primary_fold_assignment.csv", index_col=0)

y_map = {"EOCRC": 1, "LOCRC": 0}
N_TOP_GENES = 2000
MIN_CPM, MIN_FRAC_DETECTED = 1.0, 0.7

libsize = counts.sum(axis=0)
cpm = counts.div(libsize, axis=1) * 1e6

all_true, all_pred_proba, all_pred_label, all_patient = [], [], [], []
per_fold_results = []

for fold_i in sorted(folds.fold.unique()):
    train_p = folds[folds.fold != fold_i].index.tolist()
    test_p = folds[folds.fold == fold_i].index.tolist()
    train_cohort = folds.loc[train_p, "cohort"]

    eocrc_train = train_cohort[train_cohort == "EOCRC"].index
    locrc_train = train_cohort[train_cohort == "LOCRC"].index

    # unsupervised expression-robustness filter -- deliberately does not use the age
    # label, see docs/decision_log.md item 3 for why
    detected_eocrc = (cpm.loc[:, eocrc_train] > MIN_CPM).mean(axis=1) >= MIN_FRAC_DETECTED
    detected_locrc = (cpm.loc[:, locrc_train] > MIN_CPM).mean(axis=1) >= MIN_FRAC_DETECTED
    robust_genes = cpm.index[detected_eocrc & detected_locrc]

    y_train = train_cohort.map(y_map).values
    y_test = folds.loc[test_p, "cohort"].map(y_map).values

    top_genes = log2cpm.loc[robust_genes, train_p].var(axis=1).nlargest(N_TOP_GENES).index

    X_train = log2cpm.loc[top_genes, train_p].T.values
    X_test = log2cpm.loc[top_genes, test_p].T.values

    scaler = StandardScaler().fit(X_train)
    X_train_s = scaler.transform(X_train)
    X_test_s = scaler.transform(X_test)

    inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    clf = LogisticRegressionCV(
        penalty="elasticnet", solver="saga",
        l1_ratios=[0.3, 0.7], Cs=5,
        cv=inner_cv, max_iter=1000, tol=1e-3, scoring="roc_auc",
        random_state=42, n_jobs=1
    )
    clf.fit(X_train_s, y_train)

    proba = clf.predict_proba(X_test_s)[:, 1]
    pred = clf.predict(X_test_s)

    fold_auc = roc_auc_score(y_test, proba) if len(set(y_test)) > 1 else np.nan
    fold_acc = accuracy_score(y_test, pred)
    per_fold_results.append(dict(fold=int(fold_i), n_test=len(test_p), auc=fold_auc, acc=fold_acc))

    all_true.extend(y_test)
    all_pred_proba.extend(proba)
    all_pred_label.extend(pred)
    all_patient.extend(test_p)

    print(f"Fold {fold_i}: n_test={len(test_p)}, AUC={fold_auc:.3f}, Acc={fold_acc:.3f}")

overall_auc = roc_auc_score(all_true, all_pred_proba)
overall_acc = accuracy_score(all_true, all_pred_label)
print(f"\n=== Representation A, discovery cohort (n=98) ===")
print(f"Overall AUC: {overall_auc:.3f} | Overall Accuracy: {overall_acc:.3f}")

pd.DataFrame(per_fold_results).to_csv("data/repA_per_fold_results.csv", index=False)
pd.DataFrame({"patient": all_patient, "true": all_true,
              "pred_proba": all_pred_proba, "pred_label": all_pred_label}) \
    .to_csv("data/repA_oof_predictions.csv", index=False)
