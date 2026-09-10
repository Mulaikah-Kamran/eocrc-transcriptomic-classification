"""Figure 2 -- Main performance comparison: discovery vs. external TCGA AUC across
the three representations. Loads from the actual results table rather than
hardcoded numbers, so the figure always matches results/tables/table2.

Reads: results/tables/table2_representation_comparison.csv
Produces: results/figures/figure2_performance_comparison.png
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

os.makedirs("results/figures", exist_ok=True)

table2 = pd.read_csv("results/tables/table2_representation_comparison.csv")
representations = ["Gene-level\n(2,000 genes)", "Hallmark pathway\n(50 sets)", "Marx 8-gene\npanel"]
discovery_auc = table2["Discovery CV AUC"].tolist()
external_auc = table2["External TCGA AUC"].tolist()
external_p = table2["External permutation p"].tolist()

x = np.arange(len(representations))
width = 0.35

fig, ax = plt.subplots(figsize=(8, 5.5))
bars1 = ax.bar(x - width / 2, discovery_auc, width, label="Discovery CV AUC (GSE213092)", color="#4C72B0")
bars2 = ax.bar(x + width / 2, external_auc, width, label="External AUC (TCGA-COAD)", color="#DD8452")

ax.axhline(0.5, color="gray", linestyle="--", linewidth=1, label="Chance (AUC=0.5)")
ax.set_ylabel("ROC-AUC")
ax.set_title("Discovery vs. External Validation Performance by Representation")
ax.set_xticks(x)
ax.set_xticklabels(representations)
ax.set_ylim(0.4, 0.75)
ax.legend(loc="upper right", fontsize=9)

for i, (b1, b2, p) in enumerate(zip(bars1, bars2, external_p)):
    ax.text(b1.get_x() + b1.get_width() / 2, b1.get_height() + 0.008, f"{discovery_auc[i]:.3f}",
            ha="center", fontsize=9)
    ax.text(b2.get_x() + b2.get_width() / 2, b2.get_height() + 0.008, f"{external_auc[i]:.3f}",
            ha="center", fontsize=9)
    ax.text(b2.get_x() + b2.get_width() / 2, 0.42, f"perm p={p}", ha="center", fontsize=8,
            style="italic", color="dimgray")

plt.tight_layout()
plt.savefig("results/figures/figure2_performance_comparison.png", dpi=200)
print("Saved results/figures/figure2_performance_comparison.png")
