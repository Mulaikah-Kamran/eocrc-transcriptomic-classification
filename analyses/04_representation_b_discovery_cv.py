"""Step 04 -- Representation B (Hallmark pathway-level), discovery-cohort cross-validation.

Reads:
    data/primary_log2cpm.csv, data/primary_fold_assignment.csv
    data/raw/hallmark_gene_sets.gmt  (MSigDB Hallmark, gene symbols -- see docs/data_provenance.md)
Produces: data/repB_per_fold_results.csv, data/repB_oof_predictions.csv
"""
import sys
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, ".")
from src.gmt_utils import parse_gmt

log2cpm = pd.read_csv("data/primary_log2cpm.csv", index_col=0)
folds = pd.read_csv("data/primary_fold_assignment.csv", index_col=0)
genesets = parse_gmt("data/raw/hallmark_gene_sets.gmt", restrict_to=log2cpm.index)

print(f"Pathway gene-set sizes after intersecting with our data: "
      f"min={min(len(v) for v in genesets.values())}, "
      f"mean={round(sum(len(v) for v in genesets.values())/len(genesets))}, "
      f"max={max(len(v) for v in genesets.values())}")

y_map = {"EOCRC": 1, "LOCRC": 0}


def compute_pathway_scores(train_patients, all_patients):
    """z-score each gene using TRAIN patients only, then average by pathway."""
    train_expr = log2cpm[train_patients]
    gene_mean = train_expr.mean(axis=1)
    gene_std = train_expr.std(axis=1).replace(0, np.nan)
    z = log2cpm[all_patients].sub(gene_mean, axis=0).div(gene_std, axis=0)
    return pd.DataFrame({name: z.loc[list(genes)].mean(axis=0) for name, genes in genesets.items()})


all_true, all_pred_proba, all_pred_label, all_patient = [], [], [], []
per_fold_results = []

for fold_i in sorted(folds.fold.unique()):
    train_p = folds[folds.fold != fold_i].index.tolist()
    test_p = folds[folds.fold == fold_i].index.tolist()

    scores = compute_pathway_scores(train_p, train_p + test_p)
    X_train = scores.loc[train_p].values
    X_test = scores.loc[test_p].values
    y_train = folds.loc[train_p, "cohort"].map(y_map).values
    y_test = folds.loc[test_p, "cohort"].map(y_map).values

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
print(f"\n=== Representation B (Hallmark pathway z-score), discovery cohort ===")
print(f"Overall AUC: {overall_auc:.3f} | Overall Accuracy: {overall_acc:.3f}")

pd.DataFrame(per_fold_results).to_csv("data/repB_per_fold_results.csv", index=False)
pd.DataFrame({"patient": all_patient, "true": all_true,
              "pred_proba": all_pred_proba, "pred_label": all_pred_label}) \
    .to_csv("data/repB_oof_predictions.csv", index=False)
