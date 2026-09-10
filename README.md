# Comparing Transcriptomic Representations for Cross-Cohort Classification of Early- and Late-Onset Colorectal Cancer

**Can gene expression tell you whether a colorectal tumor came from a younger or older patient — and does that signal hold up when tested somewhere else?**

The first answer I got was wrong. Figuring out why turned into the most valuable part of the project.

## The short version

I built a classifier separating early- from late-onset colorectal cancer patients by gene expression. It worked almost too well — which was the giveaway. The patients had been sequenced on different machines; the model was reading lab equipment, not biology. I threw that result out, found a cleaner dataset, and rebuilt the analysis.

On the clean data: a gene-level model found a real signal that held up on an independent cohort from The Cancer Genome Atlas. A pathway-based model and a small gene panel from a previously published study did not. Looking at what genes the working model relied on pointed toward the tumor's surrounding tissue — blood vessels, connective tissue, immune cells — more than the tumor cells themselves.

Full write-up: **[REPORT.md](REPORT.md)**.

## The headline result

| Representation | Discovery cohort | Independent test (TCGA-COAD) | Better than chance? |
|---|---|---|---|
| Gene-level (2,000 genes) | AUC 0.666 | **AUC 0.626** | **Yes — p = 0.001** |
| Pathway-level (50 Hallmark gene sets) | AUC 0.624 | AUC 0.548 | No |
| Published 8-gene panel | AUC 0.591 | AUC 0.540 | No |

AUC 0.626 is modest — not a diagnostic tool. What makes it real: it reproduced in a completely independent group of patients, tested against a bar most exploratory transcriptomics studies never attempt.

## Where this came from

This builds on an earlier internship analyzing tumor-versus-normal RNA-seq in colorectal cancer — where I first got comfortable with the full pipeline, from raw sequencing data to biological interpretation. This project asks a harder question with the same disease.

## What's in this repository

- **[REPORT.md](REPORT.md)** — full write-up: background, methods, results, discussion, references.
- **`docs/`** — the decisions behind the science:
  - `data_provenance.md` — data sources and QC steps
  - `decision_log.md` — every major methodological choice and why
  - `limitations.md` — what this project does and doesn't establish
  - `literature_review.md` — the published research this is grounded in
- **`analyses/`** — the code, numbered in run order
- **`results/`** — figures and tables

## Reproducing this

```bash
pip install -r requirements.txt
```

Place these raw files in `data/raw/` (see `docs/data_provenance.md` for accessions and sources — not included here due to size and redistribution terms):

- `GSE213092_series_matrix`, `GSE213092_HRCRC.txt` (discovery cohort, GEO)
- `TCGA-COAD.clinical.tsv.gz`, `TCGA-COAD.star_counts.tsv.gz` (external validation, UCSC Xena)
- `hallmark_gene_sets.gmt` (MSigDB Hallmark, gene symbols)
- `gene_annotation_reference.csv` (Ensembl-to-symbol lookup for TCGA harmonization)

Run `analyses/` in numbered order (01–18). Note: `06_permutation_test_discovery_1200.py` runs all 1,200 permutations in one pass by default (~15 minutes) — see the script header to chunk it instead.
