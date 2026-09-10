"""Figure 3 -- Representation A interpretation: top SHAP features + bootstrap stability.

Reads: data/repA_shap_importance.csv, data/repA_gene_stability.csv
Produces: results/figures/figure3_shap_stability.png
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

os.makedirs("results/figures", exist_ok=True)

shap_imp = pd.read_csv("data/repA_shap_importance.csv", index_col=0).iloc[:, 0].sort_values(ascending=False)
stability = pd.read_csv("data/repA_gene_stability.csv", index_col=0)

top20 = shap_imp.head(20)

fig, axes = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={'width_ratios': [1.3, 1]})

axes[0].barh(top20.index[::-1], top20.values[::-1], color="#4C72B0")
axes[0].set_xlabel("Mean |SHAP value|")
axes[0].set_title("Top 20 model-associated genes\n(Representation A, frozen model)")

n_stable = (stability["nonzero_rate"] == 1.0).sum()
axes[1].hist(stability["nonzero_rate"], bins=20, color="#55A868", edgecolor="white")
axes[1].set_xlabel("Fraction of 200 bootstraps with non-zero coefficient")
axes[1].set_ylabel("Number of genes (of 2,000)")
axes[1].set_title(f"Bootstrap stability distribution\n({n_stable} genes stable in all 200 resamples)")
axes[1].axvline(1.0, color="firebrick", linestyle="--", linewidth=1)

plt.tight_layout()
plt.savefig("results/figures/figure3_shap_stability.png", dpi=200)
print("Saved results/figures/figure3_shap_stability.png")
