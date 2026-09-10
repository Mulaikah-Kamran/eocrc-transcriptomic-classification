"""Step 16 -- Hypergeometric Hallmark pathway enrichment of the high-confidence gene set.

Background = 12,865 robust-expression-filtered candidate genes (the genes that were
ever eligible to be selected by the model) -- NOT all 33,115 raw genes.
Foreground = genes with non-zero coefficient in all 200 bootstrap resamples.
FDR (Benjamini-Hochberg) applied across all 50 Hallmark pathways.

Reads:
    data/primary_counts_raw.csv, data/primary_fold_assignment.csv
    data/repA_high_confidence_genes.csv
    data/raw/hallmark_gene_sets.gmt
Produces:
    results/tables/table4_full_hallmark_enrichment.csv     (all 50 pathways)
    data/repA_pathway_enrichment_significant_with_genes.csv (FDR<0.05 subset, with gene lists)
"""
import sys
import os
import pandas as pd
import numpy as np
from scipy.stats import hypergeom

sys.path.insert(0, ".")
from src.gmt_utils import parse_gmt

os.makedirs("results/tables", exist_ok=True)

counts = pd.read_csv("data/primary_counts_raw.csv", index_col=0)
folds = pd.read_csv("data/primary_fold_assignment.csv", index_col=0)
libsize = counts.sum(axis=0)
cpm = counts.div(libsize, axis=1) * 1e6
eocrc_p = folds[folds.cohort == "EOCRC"].index.tolist()
locrc_p = folds[folds.cohort == "LOCRC"].index.tolist()
detected_e = (cpm[eocrc_p] > 1.0).mean(axis=1) >= 0.7
detected_l = (cpm[locrc_p] > 1.0).mean(axis=1) >= 0.7
background = set(cpm.index[detected_e & detected_l])

high_conf = pd.read_csv("data/repA_high_confidence_genes.csv", index_col=0)
foreground = set(high_conf.index) & background
print(f"Background: {len(background)}, Foreground: {len(foreground)}")

genesets = parse_gmt("data/raw/hallmark_gene_sets.gmt", restrict_to=background)

N, n = len(background), len(foreground)
results = []
for name, pathway_genes in genesets.items():
    K = len(pathway_genes)
    if K == 0:
        continue
    k = len(pathway_genes & foreground)
    pval = hypergeom.sf(k - 1, N, K, n)
    fold_enrichment = (k / n) / (K / N) if K > 0 and n > 0 else np.nan
    results.append(dict(pathway=name, k_in_foreground=k, K_pathway_size=K,
                         fold_enrichment=fold_enrichment, pval=pval))

res_df = pd.DataFrame(results).sort_values("pval").reset_index(drop=True)
m = len(res_df)
res_df["rank"] = np.arange(1, m + 1)
res_df["fdr_raw"] = res_df["pval"] * m / res_df["rank"]
res_df["fdr"] = res_df["fdr_raw"][::-1].cummin()[::-1].clip(upper=1.0)

print(f"Pathways tested: {m}, significant at FDR<0.05: {(res_df.fdr < 0.05).sum()}")
res_df.to_csv("results/tables/table4_full_hallmark_enrichment.csv", index=False)

sig = res_df[res_df.fdr < 0.05].copy()
sig["contributing_genes"] = sig["pathway"].apply(
    lambda name: ", ".join(sorted(genesets[name] & foreground))
)
sig.to_csv("data/repA_pathway_enrichment_significant_with_genes.csv", index=False)
print("Saved results/tables/table4_full_hallmark_enrichment.csv and "
      "data/repA_pathway_enrichment_significant_with_genes.csv")
