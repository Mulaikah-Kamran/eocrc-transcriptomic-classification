"""Figure 4 -- Biological audit: top enriched Hallmark pathways + cell-type marker
membership. Marker presence/absence is loaded from the actual analysis output
(results/tables/cell_composition_results.csv) rather than hardcoded -- an earlier
version of this script hardcoded these values directly, which meant the figure
could silently drift out of sync with the real result.

Reads:
    data/repA_pathway_enrichment_significant_with_genes.csv
    results/tables/cell_composition_results.csv
Produces: results/figures/figure4_biological_audit.png
"""
import os
import ast
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

os.makedirs("results/figures", exist_ok=True)

enrich = pd.read_csv("data/repA_pathway_enrichment_significant_with_genes.csv").sort_values("pval")
top10 = enrich.head(10).iloc[::-1]

composition = pd.read_csv("results/tables/cell_composition_results.csv")


def has_any_marker(cell_str):
    try:
        parsed = ast.literal_eval(cell_str)
    except (ValueError, SyntaxError):
        parsed = cell_str
    return len(parsed) > 0 if isinstance(parsed, list) else bool(parsed)


composition["present"] = composition["markers_present_in_high_confidence_set"].apply(has_any_marker)

fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [1.4, 1]})

colors = -np.log10(top10["fdr"])
axes[0].barh(top10["pathway"].str.replace("HALLMARK_", "").str.replace("_", " "),
             top10["fold_enrichment"], color=plt.cm.Reds(colors / colors.max()))
axes[0].set_xlabel("Fold enrichment")
axes[0].set_title("Top 10 enriched Hallmark pathways\n(high-confidence gene set vs. background)")
for i, (fe, fdr) in enumerate(zip(top10["fold_enrichment"], top10["fdr"])):
    axes[0].text(fe + 0.1, i, f"FDR={fdr:.1e}", va="center", fontsize=8)

names = composition["compartment"].tolist()
vals = composition["present"].astype(int).tolist()
bar_colors = ["#55A868" if v == 1 else "#C44E52" for v in vals]
axes[1].barh(names, [1] * len(names), color=bar_colors)
axes[1].set_xlim(0, 1.3)
axes[1].set_xticks([])
axes[1].set_title("Canonical marker presence in\nhigh-confidence gene set")
for i, v in enumerate(vals):
    axes[1].text(1.03, i, "Present" if v == 1 else "Absent", va="center", fontsize=9,
                 color="#55A868" if v == 1 else "#C44E52", fontweight="bold")

plt.tight_layout()
plt.savefig("results/figures/figure4_biological_audit.png", dpi=200)
print("Saved results/figures/figure4_biological_audit.png")
