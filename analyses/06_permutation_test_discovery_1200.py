"""Step 06 -- Extended discovery-cohort permutation test (1,200 shuffles total).

Running with no arguments does all 1,200 shuffles in one process (this will take
a while -- roughly 0.7s per shuffle on a single core, so ~15 minutes total).

To split the work across multiple runs (e.g. if you want to parallelize manually,
or resume after an interruption), pass three arguments:

    python3 analyses/06_permutation_test_discovery_1200.py <chunk_id> <n_shuffles> <seed_offset>

    Example: run 6 chunks of 200 shuffles each, on 6 different seeds:
        python3 analyses/06_permutation_test_discovery_1200.py 1 200 1
        python3 analyses/06_permutation_test_discovery_1200.py 2 200 2
        ... etc, then combine the resulting data/perm1200_chunk*.csv files.

Reads: data/primary_counts_raw.csv, data/primary_log2cpm.csv, data/primary_fold_assignment.csv
Produces: data/perm1200_chunk<chunk_id>.csv (combine all chunks for the full result)
"""
import sys
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings("ignore")

# optional CLI args, with defaults that make `python3 06_....py` runnable as-is
chunk_id = int(sys.argv[1]) if len(sys.argv) > 1 else 1
n_this_chunk = int(sys.argv[2]) if len(sys.argv) > 2 else 1200
seed_offset = int(sys.argv[3]) if len(sys.argv) > 3 else 1

log2cpm = pd.read_csv("data/primary_log2cpm.csv", index_col=0)
folds = pd.read_csv("data/primary_fold_assignment.csv", index_col=0)
counts = pd.read_csv("data/primary_counts_raw.csv", index_col=0)
libsize = counts.sum(axis=0)
cpm = counts.div(libsize, axis=1) * 1e6
y_map = {"EOCRC": 1, "LOCRC": 0}
N_TOP_GENES = 2000
true_labels = folds.cohort.map(y_map)


def run_cv(label_series):
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
        clf = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=1.0,
                                  max_iter=500, tol=1e-2, random_state=0)
        clf.fit(scaler.transform(X_train), y_train)
        proba = clf.predict_proba(scaler.transform(X_test))[:, 1]
        all_true.extend(y_test)
        all_proba.extend(proba)
    return roc_auc_score(all_true, all_proba)


rng = np.random.default_rng(1000 + seed_offset)
results = []
for i in range(n_this_chunk):
    shuffled = pd.Series(rng.permutation(true_labels.values), index=true_labels.index)
    results.append(run_cv(shuffled))
    if (i + 1) % 100 == 0:
        print(f"  {i+1}/{n_this_chunk} permutations done")

out_path = f"data/perm1200_chunk{chunk_id}.csv"
pd.Series(results, name="perm_auc").to_csv(out_path, index=False)
print(f"Chunk {chunk_id}: {n_this_chunk} permutations -> saved to {out_path}")

print("\nIf this was the only/last chunk, combine and summarize with:")
print("  python3 -c \"import pandas as pd, glob, numpy as np; "
      "aucs = pd.concat([pd.read_csv(f) for f in glob.glob('data/perm1200_chunk*.csv')]); "
      "print(len(aucs), 'total permutations, mean', aucs.perm_auc.mean())\"")
