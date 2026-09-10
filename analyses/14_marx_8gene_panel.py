"""Step 14 -- Externally nominated Marx et al. 8-gene panel: train on GSE213092,
freeze, apply to TCGA-COAD. See docs/decision_log.md for why ridge (not Elastic
Net) is used here, and why the gene list is never modified based on performance.

Reads:
    data/primary_log2cpm.csv, data/primary_fold_assignment.csv
    tcga/tcga_log2cpm_symbols.csv, tcga/tumor_final_audit.csv
Produces:
    data/frozen_marx8_model.pkl
    data/tcga_prediction_scores_marx8.csv
"""
import pandas as pd
import numpy as np
import pickle
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, confusion_matrix, balanced_accuracy_score, recall_score
import warnings
warnings.filterwarnings("ignore")

MARX_GENES = ["ALDOB", "FBXL16", "IL1RN", "MSLN", "RAC3", "SLC38A11", "WBSCR27", "WNT11"]

log2cpm = pd.read_csv("data/primary_log2cpm.csv", index_col=0)
folds = pd.read_csv("data/primary_fold_assignment.csv", index_col=0)
y_map = {"EOCRC": 1, "LOCRC": 0}

# --- discovery-cohort cross-validation ---
all_true, all_proba = [], []
for fold_i in sorted(folds.fold.unique()):
    train_p = folds[folds.fold != fold_i].index.tolist()
    test_p = folds[folds.fold == fold_i].index.tolist()
    y_train = folds.loc[train_p, "cohort"].map(y_map).values
    y_test = folds.loc[test_p, "cohort"].map(y_map).values
    X_train = log2cpm.loc[MARX_GENES, train_p].T.values
    X_test = log2cpm.loc[MARX_GENES, test_p].T.values
    scaler = StandardScaler().fit(X_train)
    clf = LogisticRegressionCV(penalty="l2", solver="lbfgs", Cs=10,
                                cv=StratifiedKFold(5, shuffle=True, random_state=42),
                                max_iter=2000, scoring="roc_auc", random_state=42)
    clf.fit(scaler.transform(X_train), y_train)
    proba = clf.predict_proba(scaler.transform(X_test))[:, 1]
    all_true.extend(y_test)
    all_proba.extend(proba)
disc_auc = roc_auc_score(all_true, all_proba)
print(f"Marx 8-gene panel -- discovery cohort CV AUC: {disc_auc:.4f}")

# --- freeze final model on all 98 patients ---
all_patients = folds.index.tolist()
y_all = folds["cohort"].map(y_map).values
X_all = log2cpm.loc[MARX_GENES, all_patients].T.values
scaler_final = StandardScaler().fit(X_all)
final_model = LogisticRegressionCV(penalty="l2", solver="lbfgs", Cs=10,
                                    cv=StratifiedKFold(5, shuffle=True, random_state=42),
                                    max_iter=2000, scoring="roc_auc", random_state=42)
final_model.fit(scaler_final.transform(X_all), y_all)

with open("data/frozen_marx8_model.pkl", "wb") as f:
    pickle.dump({"genes": MARX_GENES, "scaler_mean": scaler_final.mean_,
                 "scaler_scale": scaler_final.scale_, "model": final_model,
                 "discovery_cv_auc": disc_auc}, f)

# --- apply to TCGA (reusing the shared tcga_log2cpm_symbols.csv) ---
log2cpm_tcga = pd.read_csv("tcga/tcga_log2cpm_symbols.csv", index_col=0)
tumor_valid = pd.read_csv("tcga/tumor_final_audit.csv")

scaler_mean = pd.Series(scaler_final.mean_, index=MARX_GENES)
scaler_scale = pd.Series(scaler_final.scale_, index=MARX_GENES)
X_tcga = pd.DataFrame(index=tumor_valid["sample"], columns=MARX_GENES, dtype=float)
for g in MARX_GENES:
    if g in log2cpm_tcga.index:
        X_tcga[g] = log2cpm_tcga.loc[g, tumor_valid["sample"]].values
missing = [g for g in MARX_GENES if g not in log2cpm_tcga.index]
print(f"Marx genes missing from TCGA: {missing}")  # WBSCR27 is expected here -- see docs/limitations.md

X_tcga_z = (X_tcga - scaler_mean) / scaler_scale
X_tcga_z[missing] = 0.0

y_true = tumor_valid.set_index("sample").loc[X_tcga_z.index, "age_group"].map(y_map).values
proba = final_model.predict_proba(X_tcga_z[MARX_GENES].values)[:, 1]
pred = (proba >= 0.5).astype(int)

auc = roc_auc_score(y_true, proba)
bal_acc = balanced_accuracy_score(y_true, pred)
print(f"\n=== Marx 8-gene panel: TCGA-COAD external validation ===")
print(f"ROC-AUC: {auc:.3f} | Balanced accuracy: {bal_acc:.3f}")
print(f"Confusion matrix:\n{confusion_matrix(y_true, pred)}")

pd.DataFrame({"proba": proba, "true": np.where(y_true == 1, "EOCRC", "LOCRC")}) \
    .to_csv("data/tcga_prediction_scores_marx8.csv", index=False)
