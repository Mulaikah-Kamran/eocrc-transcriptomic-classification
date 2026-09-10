"""Step 01 -- Parse the raw GSE213092 series matrix and counts file.

Expects, per docs/data_provenance.md:
    data/raw/GSE213092_series_matrix
    data/raw/GSE213092_HRCRC.txt

Produces:
    data/gse213092_metadata_full.csv        (all 99 samples, patient-level metadata)
    data/gse213092_counts_full.csv          (all 99 samples x cleaned gene set)

This step only cleans the deposited file (removes malformed rows) and parses
metadata -- it does NOT remove the library-size outlier or build folds yet,
that happens in 02_build_primary_cohort_and_folds.py.
"""
import os
import csv
import pandas as pd

RAW_DIR = "data/raw"
OUT_DIR = "data"
os.makedirs(OUT_DIR, exist_ok=True)


def parse_series_matrix(path):
    rows = {}
    with open(path) as f:
        for line in f:
            if line.startswith("!Sample_"):
                parts = next(csv.reader([line], delimiter="\t"))
                rows.setdefault(parts[0], []).append(parts[1:])
    return rows


def clean(s):
    return s.strip('"')


sm = parse_series_matrix(os.path.join(RAW_DIR, "GSE213092_series_matrix"))
titles = [clean(x) for x in sm["!Sample_title"][0]]
gsm = [clean(x) for x in sm["!Sample_geo_accession"][0]]
tissue_site = [clean(x).split(": ", 1)[1] for x in sm["!Sample_characteristics_ch1"][0]]
age = [int(clean(x).split(": ", 1)[1]) for x in sm["!Sample_characteristics_ch1"][1]]
sex = [clean(x).split(": ", 1)[1] for x in sm["!Sample_characteristics_ch1"][2]]
disease = [clean(x).split(": ", 1)[1] for x in sm["!Sample_characteristics_ch1"][3]]

meta = pd.DataFrame({"hrcrc_id": titles, "gsm": gsm, "tissue_site": tissue_site,
                      "age": age, "sex": sex, "disease_state": disease})
meta["cohort"] = meta.disease_state.map({"Early-onset colorectal cancer": "EOCRC",
                                          "Late-onset colorectal cancer": "LOCRC"})
meta.to_csv(os.path.join(OUT_DIR, "gse213092_metadata_full.csv"), index=False)
print(f"Parsed metadata for {len(meta)} samples "
      f"({(meta.cohort=='EOCRC').sum()} EOCRC, {(meta.cohort=='LOCRC').sum()} LOCRC)")

# --- counts file: drop the 10 malformed rows present in the deposited file ---
counts = pd.read_csv(os.path.join(RAW_DIR, "GSE213092_HRCRC.txt"), sep="\t", index_col=0)
before = counts.shape[0]
counts = counts.dropna(axis=0)
counts = counts[~counts.index.duplicated(keep=False)]
print(f"Dropped {before - counts.shape[0]} malformed/duplicate rows; "
      f"{counts.shape[0]} genes remain")

counts.to_csv(os.path.join(OUT_DIR, "gse213092_counts_full.csv"))
print("Saved data/gse213092_metadata_full.csv and data/gse213092_counts_full.csv")
