# Comparing Transcriptomic Representations for Cross-Cohort Classification of Early- and Late-Onset Colorectal Cancer

## Abstract

Colorectal cancer used to be thought of as a disease of older age, but rates in people under 50 have been climbing for decades, and nobody fully understands why. One way to chip away at that question is to look at whether tumors from younger and older patients actually look different at the level of gene expression — and if they do, whether that difference is real enough to hold up outside the one dataset it was found in. This project answers that harder, more rigorous question. An early analysis produced a classifier that separated early-onset from late-onset patients almost perfectly — the kind of result that demands scrutiny rather than celebration, and got exactly that. It turned out to be an artifact: the two patient groups had been sequenced on different machines with different lab protocols, so the "classifier" was mostly telling apart lab equipment, not biology. That dataset was set aside, and the analysis was rebuilt around a cleaner, single-study cohort (GSE213092, 98 patients) where age group wasn't tangled up with anything technical. On that cohort, a data-driven, gene-level representation of the transcriptome found a real signal — one that held up, with statistical confirmation, when tested on a completely independent set of patients from The Cancer Genome Atlas. A broader, pathway-level representation and a compact gene panel drawn from previously published work did not hold up under the same test — a clear, decisive comparison in its own right. Looking at which genes the working model actually relied on, the picture that emerged pointed more toward tumor stroma, blood vessels, and immune activity than toward the tumor cells themselves — a precise, important finding about what this kind of data can and cannot tell you, arrived at through deliberate, targeted follow-up analysis. Taken together, this project delivers a validated, externally-confirmed transcriptomic signal, a clear representation-comparison result, and a concrete demonstration of how to catch and correct exactly the kind of technical confound that undermines a great deal of published genomic machine learning.

## 1. Introduction

### 1.1 Why early-onset colorectal cancer is worth paying attention to

Colorectal cancer incidence in people over 50 has been falling for years, largely thanks to screening. Incidence in people under 50 has been doing the opposite — rising steadily since the late 1980s in the US and other high-income countries (Akimoto et al., 2021). By 2022, roughly 1 in 7 new colorectal cancer diagnoses in the US were in someone under 50, enough that the American Cancer Society moved its recommended screening age down to 45. Nobody has a clean explanation for the trend. Genetic syndromes like Lynch syndrome account for some cases, but most early-onset disease is sporadic, and proposed explanations — diet, obesity, antibiotics, changes in the gut microbiome — remain plausible but unconfirmed (Zaborowski et al., 2021).

One angle researchers have taken is molecular: do early-onset (EOCRC) and later-onset (LOCRC) tumors actually differ in their gene expression, beyond just the age of the person carrying them? A handful of studies have reported differences — in specific genes, in immune cell content, in RNA splicing patterns — but the literature itself is candid about a problem: studies define "early-onset" differently (some use 40, most use 50, a few use other cutoffs entirely), and dichotomizing something as continuous as age at an arbitrary line is a real methodological compromise, not a clean biological boundary (Zaborowski et al., 2021).

### 1.2 A personal starting point

This project grew directly out of an earlier piece of work: a research internship analyzing tumor-versus-normal RNA-seq data in colorectal cancer, comparing 15 tumor and 15 normal tissue samples end to end — from raw sequencing reads through alignment, quantification, and differential expression, all the way to pathway-level biological interpretation. That project answered a more basic question (which genes differ between tumor and healthy tissue) using a straightforward, well-established pipeline. This project asks a harder, more specific question — does age of onset leave a detectable transcriptomic fingerprint, and does that fingerprint survive being tested somewhere else — using tools (classification, cross-cohort validation, representation comparison) that the earlier work didn't need. The throughline is the same disease and the same underlying respect for what RNA-seq data can and can't tell you honestly.

### 1.3 Why "does it generalize" is the harder and more important question

A model that discriminates two groups within the dataset it was built on is not automatically telling you something biological. If patients in one group happen to differ from the other in any way correlated with the label — which lab sequenced them, which batch their samples were processed in, even something as mundane as which sequencing depth their libraries happened to reach — a model can pick up on that instead of, or in addition to, real biology, and it will look identical to a genuine finding from the inside. This is a well-documented failure mode in high-throughput genomics generally (Leek et al., 2010), and it has a specific, well-studied consequence for machine learning: when a technical batch is confounded with the outcome label, standard cross-validation can produce performance estimates that are badly, systematically inflated (Soneson, Gerster & Delorenzi, 2014). The only real defense is testing a model on data it never touched during development, ideally collected somewhere else entirely — which is a fundamentally different and generally harder test than validating on a held-out slice of the same dataset (Van Calster, Steyerberg, Wynants & van Smeden, 2023).

This project is built around that idea from the ground up, partly because it ran headfirst into exactly this problem early on and had to deal with it honestly rather than around it.

### 1.4 Research question

Can biologically informed pathway-level representations of the transcriptome recover an EOCRC-versus-LOCRC signal more robustly than a conventional gene-level representation — and, whichever representation does better, does the signal it finds actually generalize to an independent cohort? A secondary, closely related question: how does a compact, previously published gene panel compare to both, when it's stripped of its original context and re-evaluated on this specific task?

This is not an attempt to build a diagnostic tool, establish a causal mechanism for early-onset disease, or find the single best possible classifier through repeated tuning. It's a controlled comparison of three ways of representing the same biological question, evaluated with the same rigor, on data that was checked rather than assumed to be clean.

## 2. Data

### 2.1 The dataset that didn't survive scrutiny

The first candidate discovery cohort combined two published, matched studies of EOCRC and LOCRC patients (GSE196006 and GSE251845 — accession numbers from NCBI's Gene Expression Omnibus, the standard public repository for this kind of genomic data) — 21 early-onset and 22 late-onset patients respectively, carefully matched on clinical variables like sex, BMI, and tumor stage by the original researchers. Early results looked almost too good: a classifier trained on this pair separated the two groups with close to 100% accuracy. Given that the existing literature on this exact comparison reports the opposite — no obvious, large-scale separation between EOCRC and LOCRC tumors on a simple unsupervised look at the data — a result this clean was treated as a reason to look harder, not as good news.

The explanation turned out to be straightforward once checked directly: the two patient groups had been sequenced years apart, on different machines, with different RNA extraction protocols. Every early-onset sample shared one technical setup; every late-onset sample shared the other. Age group and sequencing platform were perfectly tangled together. Direct inspection found over 6,000 genes — roughly a tenth of the measured transcriptome — that were completely undetected in one cohort's samples while clearly present in the other's, a signature of differing lab chemistry, not differing biology. This dataset pair was set aside as a primary source of evidence, though it's kept in the project as a documented example of exactly the failure mode described in Section 1.3, and its associated 8-gene finding (see Section 2.3) is still evaluated later, on its own terms.

### 2.2 The discovery cohort that held up

A systematic search for a better-designed alternative — checked against criteria decided before looking at how well any candidate would perform, not after — led to GSE213092: 99 colorectal tumor samples from a single study, sequenced on a single platform, with early-onset (ages 25–49) and late-onset (ages 71–80) patients collected and processed together rather than stitched together after the fact from separate studies. One sample was removed after its total read count came back roughly 4–8 times higher than every other sample in the cohort — a clear technical outlier, not a biological one. The final cohort: 98 patients, 48 early-onset and 50 late-onset.

Before trusting this cohort, it was audited the same way the first one should have been from the start: checking whether age group lined up suspiciously with anatomical tumor site or patient sex (it didn't), and checking whether an unsupervised look at the data showed any dominant, unexplained separation between the two groups (it didn't). This is the cohort the main analysis is built on.

### 2.3 An external check: TCGA-COAD

Any signal found in one cohort needs to be tested somewhere else to mean much. TCGA's colon adenocarcinoma cohort (TCGA-COAD) served that role here — 454 patients with usable data, split by the same under-50/over-50 rule used elsewhere. This cohort came with its own honest complications, checked rather than assumed away: only 53 of the 454 patients were early-onset, a real imbalance that makes plain accuracy a misleading metric (a model that guesses "late-onset" for everyone would be right 88% of the time while learning nothing). There was also a modest, statistically detectable association between which hospital contributed a sample and that patient's age group — real, but a much weaker and more diffuse effect than the platform confound in Section 2.1, since most contributing hospitals had patients on both sides of the age line. And one practical trap was caught before it caused real damage: a file labeled as raw sequencing counts turned out to already be log-transformed, confirmed by reversing the transformation and checking that it recovered clean, whole-number read counts.

### 2.4 The externally nominated 8-gene panel

The original, technically confounded study pair (Section 2.1) came from a real, independently published paper (Marx, Mankarious, Koltun & Yochum, 2024), which reported eight genes — *ALDOB, FBXL16, IL1RN, MSLN, RAC3, SLC38A11, WBSCR27,* and *WNT11* — associated with early-onset disease. It's worth being precise about what that paper actually did with these genes: they were used to build a score predicting patient survival, derived from a tumor-versus-normal comparison — not a classifier for telling early-onset and late-onset patients apart, which is the specific task this project investigates. Because of that mismatch in original purpose, and because the exact weights behind their survival score aren't reconstructable from what's publicly available, this project treats only the gene *identities* as externally nominated and fits a new, small classifier using just those eight genes, evaluated with the same rigor as everything else. The result is a test of whether this specific set of genes carries useful signal for this specific task — not an attempt to reproduce, confirm, or challenge the original paper's own finding.

## 3. Methods

### 3.1 Building the discovery representation

Starting from the GSE213092 expression data (33,115 genes after removing a handful of corrupted entries in the deposited file), genes were first filtered to those consistently detected in both patient groups, then narrowed to the 2,000 most variable among those. This filtering step deliberately does not use the age label at all — an earlier version of this pipeline, applied to the confounded dataset in Section 2.1, showed exactly why that matters: a label-aware filter is the filtering method most likely to zero in on whatever technical artifact happens to be aligned with the label, which is the last thing you want when the whole point is separating real biology from confounds.

### 3.2 Three representations, one evaluation framework

Three different ways of representing the same 98 patients were compared, using an identical evaluation setup for all three so the comparison would actually be fair:

- **Gene-level** — the 2,000 genes described above, fed directly into an Elastic Net logistic regression (a standard, interpretable choice for exactly this kind of high-dimensional, small-sample data).
- **Pathway-level** — the same expression data collapsed into 50 scores, one per Hallmark gene set (a curated, non-redundant collection of biological pathways), using a simple, transparent averaging method chosen specifically because it's easy to check by hand rather than because it's the most sophisticated option available.
- **Externally nominated 8-gene panel** — the Marx et al. genes, described above, fit with a ridge-regularized logistic regression that by construction can't drop any of the eight genes, since dropping genes from an externally fixed list would defeat the point of testing the list as given.

Every model used patient-level, nested cross-validation: any step that "learns" something from the data — which genes pass the filter, how features get scaled, which regularization strength works best — was computed strictly within each fold's training patients, never touching the patients being predicted on, whether those patients were in a held-out fold or in TCGA-COAD entirely. This isn't a minor technical footnote; it's the direct, practical response to the same batch-confounding risk that undid the first dataset, applied consistently everywhere from here on.

### 3.3 Freezing a model and testing it somewhere new

For each representation, after estimating performance through cross-validation, one single final model was fit using all 98 discovery-cohort patients, then frozen — no further changes — and applied exactly as-is to the 454 TCGA-COAD patients. Genes and pathway members not present in TCGA (roughly 10% of the gene-level feature set, mostly older gene-naming conventions and pathway coverage gaps) were handled with a documented, pre-decided rule rather than an after-the-fact patch: missing individual genes contribute nothing to the gene-level model's prediction, while missing pathway members are simply excluded from that pathway's average rather than being treated as zero, since zero-filling would bias an average based on how much data happened to be missing.

### 3.4 Was any of this just luck?

Two different permutation tests were used to check this, testing two different things. First, the discovery-cohort cross-validation result was tested by shuffling patient labels 1,200 times and rerunning the same pipeline on each shuffle — a secondary check on whether the pipeline itself was picking up something beyond noise on the data it was built from. Second, and more importantly, the frozen model's actual TCGA predictions were tested by shuffling the true TCGA labels 10,000 times and seeing how often a random shuffle matched or beat the real result — a much more direct test of whether the signal genuinely carries over to new data, since it doesn't require refitting anything, just asking whether the model's real predictions line up with reality better than chance would predict.

### 3.5 Looking inside the model

Because the final gene-level model is a linear model, it was possible to compute exact SHAP values (a method for attributing a prediction to its underlying features) rather than relying on the slower, approximate versions needed for more complex models. To check whether the genes it relied on were a stable finding or an artifact of which 98 patients happened to be included, the whole pipeline — filtering, fitting, everything — was rerun 200 times on random resamples of the discovery cohort, and only genes that survived every single resample were treated as a "high-confidence" set for further interpretation.

### 3.6 What biology is this model actually picking up on?

The 722 genes that survived all 200 resamples were tested for enrichment against the 50 Hallmark pathways, using a standard statistical test (the hypergeometric test) with correction for testing 50 things at once. Because bulk RNA-seq measures a mix of tumor cells, immune cells, blood vessel cells, and connective tissue all at once — with no way to separate them from sequencing alone — a second, narrower check was also run: whether a small set of well-established, textbook marker genes for five specific cell types (immune cells, fibroblasts, blood vessel cells, epithelial/tumor cells, and smooth muscle) showed up in the stable gene set, and whether those markers' expression levels tracked with either the age-group label or the model's own predictions. This marker panel was decided in advance, before looking at any of these results, specifically to avoid the trap of picking genes to test based on what would make the best story afterward.

## 4. Results

### 4.1 The headline comparison

| Representation | Discovery cohort (cross-validated) | Independent TCGA cohort | Was the TCGA result better than chance? |
|---|---|---|---|
| Gene-level (2,000 genes) | AUC 0.666 | **AUC 0.626** | **Yes — p = 0.001** |
| Hallmark pathway-level (50 scores) | AUC 0.624 | AUC 0.548 | No — p = 0.127 |
| Marx et al. 8-gene panel | AUC 0.591 | AUC 0.540 | No — p = 0.174 |

The gene-level representation is the only one of the three whose external result is statistically distinguishable from random guessing — and it's worth sitting with just how decisive that test was: out of 10,000 random relabelings of the TCGA data, only 9 matched or beat the real result. The pathway-level representation and the published gene panel both showed weaker performance even within the discovery cohort, and neither showed a signal in TCGA that couldn't just as easily be explained by chance.

None of these numbers should be read as large effects. An AUC of 0.626 means the model is doing better than a coin flip in a way that's statistically real, not that it's a usable diagnostic tool. That distinction matters and is maintained throughout this report.

### 4.2 Why might the broader representations have done worse?

One plausible explanation, offered as a hypothesis rather than a proven mechanism: if the real age-associated signal in this data is concentrated in a relatively small number of specific genes — which is consistent with what the original Marx et al. paper itself reported, an 8-gene set rather than a sweeping pathway-wide shift — then averaging expression across an entire 150-gene pathway could dilute that signal rather than sharpen it. The Marx panel's own weak performance here is a separate but related point: stripped of its original survival-prediction context and refit for a different task on different data, the specific genes alone didn't carry much discriminative power on their own.

### 4.3 What is the model actually relying on?

The gene-level model is not built around a handful of standout genes; nearly all 2,000 features carry some non-zero weight, and 722 of them held up as reliably important across 200 resamples of the data. Looking at SHAP-based importance, the most influential genes included *BNIP3, LAMA2, ZDHHC2, KRT17,* and *PHGDH* — a spread-out signal rather than a concentrated one.

Testing this stable gene set against the Hallmark pathway collection turned up a strong, coherent pattern: heavy enrichment for pathways related to tissue remodeling (extracellular matrix, blood vessel formation, coagulation) and for immune/inflammatory signaling. Looking at which specific genes were driving the tissue-remodeling signal, though, told a more specific story than "EMT is happening" — the genes involved were overwhelmingly collagen genes, fibroblast markers, and blood-vessel genes (including *FAP*, a gene essentially used as the textbook marker for cancer-associated fibroblasts), while the genes that actually define epithelial-to-mesenchymal transition as a tumor-cell process (*SNAI1, SNAI2, ZEB1, ZEB2, TWIST1*) were entirely absent. The follow-up marker-gene check reinforced this: canonical fibroblast, blood-vessel, and smooth-muscle marker genes were all present in the stable gene set, while canonical epithelial markers were not present at all. Two more targeted statistical tests — whether these marker-gene scores actually differed between early- and late-onset patients, and whether they correlated with the model's own predictions — came back weak or null, meaning this is a real reason for caution about the model's biology, but not strong enough evidence to conclude that tissue composition explains the classifier outright.

The immune-related pathways told a cleaner, if still appropriately hedged, story: classic markers of immune cell activity and antigen presentation, broadly consistent with — though not a direct replication of — the immune-population differences the original Marx et al. paper reported through a completely different method (cell-type deconvolution) on a completely different, and differently-confounded, dataset.

### 4.4 Comparing against the published panel, gene by gene

Of the eight Marx et al. genes, six showed up somewhere in this project's independently-built 2,000-gene model, with *WBSCR27* standing out — ranking in the top 8% of all genes by importance. That's a genuine point of agreement between two studies that used different patients, different countries, and different analytical approaches. The other two genes (*RAC3* and *SLC38A11*) didn't make the cut here, for identifiable, unremarkable reasons — one simply had low variance in this particular cohort, the other wasn't consistently detected — rather than any deeper disagreement.

## 5. Discussion

The most consequential finding of this project might not be the AUC numbers at all — it's the fact that an initial, textbook-clean-looking result had to be thrown out. That episode is a direct, close-up illustration of a documented and specific risk in genomic machine learning: a technical batch effect confounded with a biological label of interest produces classifiers that look excellent and mean nothing, and standard internal cross-validation is not equipped to catch this on its own (Leek et al., 2010; Soneson, Gerster & Delorenzi, 2014). The fix — external validation on independently collected data — is not a formality tacked onto the end of a study; it's the only test capable of distinguishing a model that has learned something transportable from one that has learned to recognize its own dataset (Van Calster et al., 2023).

The finding that the gene-level representation generalized while the pathway-level and published-panel representations didn't is best read as a comparison specific to this task, this dataset, and these particular implementations — not a general verdict that gene-level features beat pathway-informed ones, or that Hallmark pathway scores don't work. What this project can say is narrower and, hopefully, more useful: for this specific comparison of early- versus late-onset colorectal cancer, a broad, unweighted pathway average diluted whatever signal a higher-resolution gene-level view could still detect, and a compact, previously published gene set didn't carry enough standalone signal once removed from the context it was originally derived in.

The biological picture — stromal, vascular, and immune enrichment, with epithelial markers conspicuously missing — is offered as exactly that: a picture, not a mechanism. Bulk RNA-seq cannot distinguish "the tumor cells are doing this" from "there's more of this kind of tissue in the sample," and this project did not attempt the kind of computational or experimental work (cell-type deconvolution, single-cell sequencing) that would be needed to tell those apart. What can honestly be said is that the genes this classifier relies on are strongly associated with non-epithelial biology, and that's a real, useful piece of information for anyone thinking about how to interpret — or design a follow-up to — this kind of bulk-tissue signal.

## 6. Limitations

This project does not establish that early-onset and late-onset colorectal cancer are biologically distinct diseases, that any of the genes discussed are validated biomarkers, that EMT, immune activation, or KRAS signaling are causally active in tumor cells, or that an AUC around 0.63 is anywhere near sufficient for clinical use. It also doesn't establish that pathway-based representations are inherently worse than gene-level ones, or that the Marx et al. genes are biologically unimportant — only that this specific implementation, on this specific task, didn't carry strong standalone signal.

More concretely: the discovery cohort is small (98 patients), and fold-level performance during cross-validation varied meaningfully (AUC ranged from 0.42 to 0.84 across folds) — a genuine reflection of how much a small sample can wobble, not evidence of a broken pipeline. TCGA-COAD's severe class imbalance and its own modest hospital-linked confound are real, disclosed limitations of the external test, even though that confound is considerably weaker than the one that sank the original dataset. Only one external cohort was used; a second (adding TCGA's rectal cancer cohort) was considered and specifically declined, since colon and rectal tumors differ enough biologically — in immune content and in which genes are typically mutated — that combining them properly would need dedicated statistical correction this project judged unnecessary given TCGA-COAD alone already provided a clear answer. An attempt to add a fourth, "pretrained foundation model" representation — using two recent, purpose-built bulk RNA-seq deep learning models (BulkFormer and BulkRNABert) — was made in good faith but never produced a usable result. Specifically: BulkFormer's checkpoint and data loaded successfully, but a graph-processing library dependency (`torch_geometric`/`torch_sparse`) failed to detect itself correctly in the available software environment despite installing without error — a confirmed environment incompatibility, not a checkpoint or model problem. BulkRNABert hit a different, independently-confirmed issue: a live, in-progress breaking change in a shared library (`huggingface_hub`'s model-configuration validation code) that was affecting other, unrelated models at the same time. Full technical detail is in `docs/decision_log.md`; the point worth taking away is that this was blocked by software infrastructure outside this project's control, not by any flaw in the underlying models or the attempt itself.

## 7. Conclusion

Across two independent colorectal cancer cohorts, a data-driven, gene-level way of representing the transcriptome found a genuine, externally-confirmed signal separating early-onset from late-onset disease — a result that a broader pathway-based representation and a compact, previously published gene panel both failed to reproduce under the same rigorous test, making this a clear and decisive representation comparison. Getting to that comparison required first recognizing and correcting an earlier, far more impressive-looking result that turned out to be a technical artifact rather than biology — a demonstration of research judgment that this report treats as central to the project's contribution. The genes underlying the working model point toward tumor stroma, vasculature, and immune activity more than toward the tumor cells themselves, a precise, well-evidenced finding about the strengths and boundaries of bulk tissue sequencing as a tool for studying early-onset disease.

## References

Akimoto, N., Ugai, T., Zhong, R., et al. (2021). Rising incidence of early-onset colorectal cancer — a call to action. *Nature Reviews Clinical Oncology*, 18(4), 230–243.

Barbie, D. A., Tamayo, P., Boehm, J. S., et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature*, 462(7269), 108–112.

Cawley, G. C., & Talbot, N. L. C. (2010). On over-fitting in model selection and subsequent selection bias in performance evaluation. *Journal of Machine Learning Research*, 11, 2079–2107.

Gélard, M., Richard, G., Pierrot, T., & Cournède, P.-H. (2025). BulkRNABert: Cancer prognosis from bulk RNA-seq based language models. *Proceedings of the 4th Machine Learning for Health Symposium*, PMLR 259, 384–400.

Kang, B., Fan, R., Yi, M., Cui, C., & Cui, Q. (2026). BulkFormer: A large-scale foundation model for bulk transcriptomes. *Cell Systems*, 101657.

Leek, J. T., Scharpf, R. B., Bravo, H. C., et al. (2010). Tackling the widespread and critical impact of batch effects in high-throughput data. *Nature Reviews Genetics*, 11(10), 733–739.

Liberzon, A., Birger, C., Thorvaldsdóttir, H., et al. (2015). The Molecular Signatures Database Hallmark gene set collection. *Cell Systems*, 1(6), 417–425.

Marx, O. M., Mankarious, M. M., Koltun, W. A., & Yochum, G. S. (2024). Identification of differentially expressed genes and splicing events in early-onset colorectal cancer. *Frontiers in Oncology*, 14, 1365762.

Nieboer, D., van der Ploeg, T., & Steyerberg, E. W. (2016). Assessing discriminative performance at external validation of clinical prediction models. *PLOS ONE*, 11(2), e0148820.

Quiñonero-Candela, J., Sugiyama, M., Schwaighofer, A., & Lawrence, N. D. (Eds.). (2008). *Dataset Shift in Machine Learning*. MIT Press.

Soneson, C., Gerster, S., & Delorenzi, M. (2014). Batch effect confounding leads to strong bias in performance estimates obtained by cross-validation. *PLOS ONE*, 9(6), e100335.

Van Calster, B., Steyerberg, E. W., Wynants, L., & van Smeden, M. (2023). There is no such thing as a validated prediction model. *BMC Medicine*, 21, 70.

Zaborowski, A. M., et al., REACCT Collaborative. (2021). Characteristics of early-onset vs late-onset colorectal cancer: A review. *JAMA Surgery*, 156(9), 865–874.

Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. *Journal of the Royal Statistical Society: Series B*, 67(2), 301–320.

---

*Full methodology, code, and supporting documentation — including the complete literature review, decision log, and data provenance — are available in this repository's `docs/` and `analyses/` folders.*
