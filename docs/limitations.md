# Limitations — What This Study Does and Does Not Establish

Precisely scoping a study's claims is a sign of rigor, not weakness — it's exactly what separates a defensible research finding from an overstated one. This section states plainly what this project's evidence supports and what it doesn't, so every claim made elsewhere in this repository can be trusted at face value.

## Does NOT establish
- Causal EMT activation in tumor cells
- That the 722 stable genes are validated biomarkers
- That KRAS signaling is activated in malignant cells
- That immune infiltration causes the classification signal
- That cell-type composition fully explains the classifier's output
- That pathway-level representations are inherently inferior to gene-level ones (only that this specific broad Hallmark mean-score representation underperformed on this specific task)
- That the Marx et al. genes are biologically irrelevant (only that this specific gene set, refit for a different endpoint, did not transfer strongly)
- That the model is clinically predictive or ready for any diagnostic use
- That AUC 0.626 is sufficient for any practical application
- That the identified genes are causal drivers of EOCRC

## DOES establish, within the stated limits
- A gene-level representation trained on a verified, non-confounded discovery cohort retains modest, statistically detectable discrimination in an independent external cohort
- Broad Hallmark pathway averaging and an externally nominated 8-gene panel did not retain comparable external discrimination under the same evaluation framework
- An initial ~100% classification result was technically confounded (cohort aliased with sequencing platform) and was correctly identified and rejected before being reported as a finding
- The genes underlying the working model are significantly enriched for stromal/vascular and immune/inflammatory biological programs, which is a real, quantified property of the model-associated gene set
- Canonical epithelial markers are notably absent from, and canonical stromal/vascular/muscle markers are notably present in, the stable gene set — a concrete reason for caution about tumor-cell-intrinsic interpretation

## Specific technical limitations

- **Discovery cohort size (n=98)** is small by machine-learning standards; fold-level AUC varied 0.42–0.84 across the 5 discovery folds, reflecting genuine sampling variability at this sample size, not instability in the pipeline.
- **TCGA-COAD class imbalance** (53 EOCRC / 401 LOCRC) means accuracy is not an informative metric; ROC-AUC was used as primary throughout.
- **A partial site confound in TCGA-COAD** (site vs. age-group association, p=0.0045 across 24 hospitals) was detected and disclosed — real, but not the complete aliasing that invalidated the originally rejected dataset.
- **Only one external cohort** was used; TCGA-READ was evaluated and explicitly not added (see decision log) due to documented biological heterogeneity between colon and rectal tumors.
- **Bulk RNA-seq cannot separate tissue composition from tumor-cell-intrinsic signaling** — the biological interpretation throughout is stated as "associated with" cell-type marker programs, never as measured cell proportions or confirmed causal biology.
- **The Hallmark mean-z-score representation is one specific pathway-scoring method**, not a definitive test of "pathway-informed ML" in general; ssGSEA or other scoring methods were not evaluated (documented as a deliberate scope decision, not an oversight).
- **Foundation-model representation (BulkFormer, BulkRNABert) evaluation was blocked by independently-verified, currently-active third-party software infrastructure issues**, not by model invalidity — no quantitative Representation C result exists.
- **WBSCR27** (one of the 8 Marx-panel genes, and the one gene individually most associated with the gene-level model) was absent from the TCGA-mapped gene space; its contribution to the Marx-panel TCGA evaluation was fixed at zero by documented convention.
- **The enrichment and cell-composition analyses are exploratory functional characterization**, not confirmatory hypothesis tests — the 722-gene foreground was itself derived from the classifier, and the 20 FDR-significant Hallmark pathways are not 20 independent findings (they substantially overlap in gene membership).
- All data are public, retrospective, and observational; no new wet-lab or clinical validation was performed.
