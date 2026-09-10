# Decision Log

Each entry: the decision, the reasoning, and what was explicitly rejected or deferred.

**1. Rejected the GSE251845+GSE196006 pair as primary discovery cohort.**
A ~100% classification result was investigated rather than accepted. Found a complete confound between age group and sequencing platform/protocol. Rejected the dataset rather than the finding — the dataset made any biological conclusion uninterpretable.

**2. Adopted GSE213092 as the primary discovery cohort.**
Selected after a systematic search against pre-specified, performance-independent criteria (direct EOCRC/LOCRC labels, genome-wide bulk RNA-seq, adequate size, single-platform consistency, no severe confounding, public accessibility, sufficient gene overlap with TCGA). Verified clean via its own confound audit (no significant tissue-site or sex association with age group; no dominant unsupervised technical separation) before adoption.

**3. Used unsupervised variance-based gene selection, not differential-expression-based selection, for Representation A.**
Chosen specifically because a supervised (label-using) filter would have been the filtering method most likely to concentrate a platform/cohort confound into the selected features — a risk directly informed by what went wrong with the rejected dataset. Kept as the standing choice for GSE213092 even after that specific risk no longer applied, for methodological consistency across the sprint.

**4. Used mean Hallmark gene-set z-scores, not ssGSEA, for Representation B.**
Chosen for transparency and hand-verifiability over a more complex rank-based algorithm that could not be independently validated against a reference implementation in this environment. Documented as an implementation choice, not a claim that this is the only or best pathway-scoring method.

**5. Patient-level, nested cross-validation throughout.**
All feature selection, normalization baselines, and hyperparameter tuning computed strictly within training folds. No patient's data ever appears in both training and test/external sets.

**6. TCGA-COAD adopted as external validation; TCGA-READ rejected.**
COAD-alone judged sufficient given it already provided a statistically informative independent result; combining with READ would require formal batch correction (standard in the literature for this specific combination) for an unclear scientific gain — assessed on non-performance criteria, not because it might change the AUC.

**7. Ridge (L2), not Elastic Net, for the Marx 8-gene panel.**
Elastic Net's L1 component could zero out some of the 8 externally-mandated genes entirely, which would silently violate the constraint of not removing genes from an externally nominated panel.

**8. Missing-gene handling differs by representation, deliberately.**
Representation A (independent gene features): missing gene → standardized value fixed at 0 (contributes nothing). Representation B (pathway averages): missing gene → excluded from that pathway's average entirely, since zero-imputing into an average would bias the score based on how many genes happened to be missing. Documented as representation-specific modeling conventions, not a single blanket rule.

**9. Foundation-model representation (BulkFormer, BulkRNABert) investigated but excluded from quantitative comparison.**
Both blocked by independently verified, currently-active third-party infrastructure incompatibilities (a `torch_geometric`/`torch_sparse` detection failure; a live `huggingface_hub` config-validation breaking change affecting unrelated models too). Not a checkpoint-validity or model-quality issue. Documented in full rather than omitted or faked with a substitute representation.

**10. Enrichment and cell-marker analyses treated as exploratory biological auditing, not confirmatory hypothesis testing.**
The 722-gene foreground was itself derived from the classifier; the 20 FDR-significant Hallmark pathways are not treated as 20 independent biological discoveries, both because of this selection dependency and because the pathways themselves overlap in gene membership (verified via pairwise Jaccard overlap before grouping into biological programs).

**11. Two final additional analyses (Marx 8-gene panel, cell-composition marker check) were pre-specified in full — including exact preprocessing, classifier choice, and what would count as an interesting result — before being run.**
Explicit guard against post-hoc rationalization of whichever result happened to look better.

**12. Hard-stopped dataset shopping after a second, targeted search found no new independent, genome-wide bulk RNA-seq EOCRC/LOCRC dataset released 2024–2026.**
Decision made independent of whether a new dataset might improve reported performance.

**13. Declined to add further classifiers, pathway databases, or a continuous-age analysis after the three-way representation comparison was complete.**
A modest AUC was treated as a legitimate scientific result requiring accurate reporting and correct framing (via literature grounding and honest limitations), not as evidence that more modeling was needed.
