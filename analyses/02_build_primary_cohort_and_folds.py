"""Step 02 -- Remove the library-size outlier and build the final 98-patient cohort,
including the locked patient-level fold assignment reused by every representation.

Reads:
    data/gse213092_metadata_full.csv
    data/gse213092_counts_full.csv

Produces:
    data/primary_metadata.csv
    data/primary_counts_raw.csv
    data/primary_log2cpm.csv
    data/primary_fold_assignment.csv
"""
import os
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold

OUTLIER_SAMPLE = "HRCRC3109"  # EOCRC sample, 186.8M reads vs. everyone else under 43M -- see docs/decision_log.md

meta = pd.read_csv("data/gse213092_metadata_full.csv")
counts = pd.read_csv("data/gse213092_counts_full.csv", index_col=0)

meta = meta[meta.hrcrc_id != OUTLIER_SAMPLE].reset_index(drop=True)
counts = counts.drop(columns=[OUTLIER_SAMPLE])
print(f"Removed outlier {OUTLIER_SAMPLE}. Final cohort: {len(meta)} patients "
      f"({(meta.cohort=='EOCRC').sum()} EOCRC, {(meta.cohort=='LOCRC').sum()} LOCRC)")

meta.to_csv("data/primary_metadata.csv", index=False)
counts.to_csv("data/primary_counts_raw.csv")

libsize = counts.sum(axis=0)
cpm = counts.div(libsize, axis=1) * 1e6
log2cpm = np.log2(cpm + 1)
log2cpm.to_csv("data/primary_log2cpm.csv")

# --- locked, patient-level, stratified fold assignment -- reused by every representation ---
labels = meta.set_index("hrcrc_id").cohort
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
fold_assignment = pd.Series(index=labels.index, dtype=int, name="fold")
for i, (_, test_idx) in enumerate(skf.split(labels.index, labels.values)):
    fold_assignment.iloc[test_idx] = i
fold_df = fold_assignment.to_frame().join(labels)
fold_df.to_csv("data/primary_fold_assignment.csv")

print("\nFold balance:")
print(fold_df.groupby(["fold", "cohort"]).size().unstack())
print("\nSaved data/primary_metadata.csv, primary_counts_raw.csv, "
      "primary_log2cpm.csv, primary_fold_assignment.csv")
