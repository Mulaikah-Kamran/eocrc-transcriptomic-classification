"""Step 11 -- Build the TCGA-COAD log2(CPM+1) matrix, symbol-indexed and harmonized.

This is the ONE place this processing happens -- scripts 12, 13, and 14 all reuse
this output rather than re-deriving it themselves (an earlier version of this
pipeline duplicated this logic in every downstream script; fixed here).

Reads:
    tcga/tumor_final_audit.csv
    data/raw/TCGA-COAD.star_counts.tsv.gz
    data/raw/gene_annotation_reference.csv
    data/harmonized_gene_renames.csv
Produces:
    tcga/tcga_log2cpm_symbols.csv
"""
import sys
import pandas as pd

sys.path.insert(0, ".")
from src.tcga_utils import reverse_xena_log2_transform, map_ensembl_to_symbol, \
    counts_to_log2cpm, apply_verified_renames

gene_info = pd.read_csv("data/raw/gene_annotation_reference.csv")
renames = pd.read_csv("data/harmonized_gene_renames.csv", index_col=0)["0"].to_dict()
tumor_valid = pd.read_csv("tcga/tumor_final_audit.csv")

log2counts_tcga = pd.read_csv("data/raw/TCGA-COAD.star_counts.tsv.gz", sep="\t", index_col=0)
keep_samples = [s for s in tumor_valid["sample"] if s in log2counts_tcga.columns]
log2counts_sub = log2counts_tcga[keep_samples].copy()

# The file is labeled "counts" but is actually log2(count+1) -- caught during the
# original audit by reversing the transform and confirming clean integer recovery.
raw_counts = reverse_xena_log2_transform(log2counts_sub)
raw_counts_by_symbol = map_ensembl_to_symbol(raw_counts, gene_info)

log2cpm_tcga = counts_to_log2cpm(raw_counts_by_symbol)
log2cpm_tcga = apply_verified_renames(log2cpm_tcga, renames)

log2cpm_tcga.to_csv("tcga/tcga_log2cpm_symbols.csv")
print(f"Saved tcga/tcga_log2cpm_symbols.csv: {log2cpm_tcga.shape[0]} genes x "
      f"{log2cpm_tcga.shape[1]} patients")
