"""Figure 1 -- Study design schematic, including the rejected dataset path.
Produces: results/figures/figure1_study_design.png
No data dependencies -- this is a fixed schematic of the study design itself.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

os.makedirs("results/figures", exist_ok=True)

fig, ax = plt.subplots(figsize=(11, 12))
ax.set_xlim(0, 10)
ax.set_ylim(0, 24)
ax.axis("off")


def box(x, y, w, h, text, color="#EAF2F8", edge="#2E75B6", fontsize=9.5, fontweight="normal"):
    b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", linewidth=1.3,
                        edgecolor=edge, facecolor=color)
    ax.add_patch(b)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize,
             fontweight=fontweight, wrap=True)


def arrow(x1, y1, x2, y2, color="black"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14,
                                  color=color, linewidth=1.2))


box(3, 22.5, 4, 1.2, "Research question:\nEOCRC vs LOCRC transcriptomic\nsignal across representations", "#D6EAF8", fontweight="bold")
arrow(5, 22.5, 5, 21.7)

box(3, 20.5, 4, 1.0, "GSE251845 + GSE196006\n(initial candidate cohort)")
arrow(5, 20.5, 5, 19.7)
box(3, 18.5, 4, 1.0, "~100% classification\n(red flag)", "#FADBD8", "#C0392B")
arrow(5, 18.5, 5, 17.7)
box(3, 16.5, 4, 1.0, "Confound audit: age group\naliased with platform/protocol", "#FADBD8", "#C0392B")
arrow(5, 16.5, 5, 15.7)
box(3, 14.5, 4, 1.0, "RESULT REJECTED\n(retained only as secondary\nhistorical comparison)", "#F5B7B1", "#943126", fontweight="bold")

arrow(5, 14.5, 5, 13.7)
box(3, 12.5, 4, 1.0, "Dataset search (performance-\nindependent criteria) -> GSE213092", "#D5F5E3", "#1E8449")
arrow(5, 12.5, 5, 11.7)
box(2.5, 10.3, 5, 1.2, "GSE213092: 98 patients\n(48 EOCRC, 50 LOCRC), QC audit\npassed (no dominant confound)")
arrow(5, 10.3, 5, 9.5)

box(0.3, 7.8, 3, 1.3, "Representation A:\n2,000 gene-level\nfeatures", "#EBF5FB")
box(3.6, 7.8, 3, 1.3, "Representation B:\n50 Hallmark pathway\nscores", "#EBF5FB")
box(6.9, 7.8, 3, 1.3, "Marx 8-gene panel\n(externally nominated)", "#EBF5FB")
arrow(1.8, 10.3, 1.8, 9.1)
arrow(5.1, 10.3, 5.1, 9.1)
arrow(8.4, 10.3, 8.4, 9.1)

for x in [1.8, 5.1, 8.4]:
    arrow(x, 7.8, x, 7.0)
box(2, 5.8, 6, 1.0, "Nested patient-level CV -> freeze final model\n(all 98 GSE213092 patients)")
arrow(5, 5.8, 5, 5.0)
box(2, 3.8, 6, 1.0, "External validation: TCGA-COAD\n(454 patients, 53 EOCRC / 401 LOCRC)")
arrow(5, 3.8, 5, 3.0)
box(2, 1.8, 6, 1.0, "Compare representations + permutation\ntesting + biological audit of stable genes", "#D6EAF8")

plt.tight_layout()
plt.savefig("results/figures/figure1_study_design.png", dpi=200)
print("Saved results/figures/figure1_study_design.png")
