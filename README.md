# Comparing Transcriptomic Representations for Cross-Cohort Classification of Early- and Late-Onset Colorectal Cancer

**Can gene expression tell you whether a colorectal tumor came from a younger or older patient — and if it can, does that signal actually hold up when you test it somewhere else?**

That question turned into a much more interesting project than expected, mostly because the first answer I got was wrong, and figuring out *why* it was wrong ended up being the most important part of the whole thing.

## The short version

I built a classifier that told apart early-onset from late-onset colorectal cancer patients using their tumors' gene expression, and it worked — almost too well, at first. It turned out that "almost too well" was the giveaway: the patients had been sequenced on different machines, so the model was mostly recognizing lab equipment, not biology. I threw that result out, found a better dataset where that couldn't happen, and rebuilt the analysis from there. On the cleaner data, a data-driven, gene-level view of the transcriptome found a real but modest signal — one that survived being tested on a completely different set of patients from The Cancer Genome Atlas. A broader, pathway-based way of representing the same data, and a small gene panel borrowed from a previously published study, both failed that same test. Looking at what genes the working model actually leaned on pointed more toward the tumor's surrounding tissue — blood vessels, connective tissue, immune cells — than toward the tumor cells themselves, which matters for how much you can actually claim from this kind of data.

The full story, written out properly, is in **[REPORT.md](REPORT.md)**.

## Where this came from

This project grew out of an earlier internship analyzing tumor-versus-normal RNA-seq data in colorectal cancer — a more straightforward comparison, but it's where I first got comfortable with the whole pipeline, from raw sequencing data through to biological interpretation. This project asks a harder question with the same disease and, I'd like to think, the same care about not overstating what the data can actually support.

## The headline result

| Representation | Discovery cohort | Independent test cohort (TCGA-COAD, a large public cancer genomics resource) | Better than chance? |
|---|---|---|---|
| Gene-level (2,000 genes) | AUC 0.666 | **AUC 0.626** | **Yes — p = 0.001** |
| Pathway-level (50 Hallmark gene sets) | AUC 0.624 | AUC 0.548 | No |
| Published 8-gene panel | AUC 0.591 | AUC 0.540 | No |

An AUC (a standard measure of how well a model separates two groups, from 0.5 = random guessing to 1.0 = perfect separation) of 0.626 isn't a diagnostic tool — it's a coin flip nudged in the right direction, and I want to be upfront about that rather than let a table of numbers imply more than it should. What makes it worth reporting is that it's a *real*, statistically checked signal that shows up again in a completely separate group of patients, which is a meaningfully higher bar than most exploratory transcriptomics analyses ever get tested against.

## What's actually in this repository

- **[REPORT.md](REPORT.md)** — the full write-up: background, methods, results, discussion, limitations, references. Written to be readable by someone without a bioinformatics background, not just other researchers.
- **`docs/`** — the paper trail behind every decision:
  - `data_provenance.md` — exactly where the data came from and every quality-control step applied to it
  - `decision_log.md` — every major methodological choice, and why it was made
  - `limitations.md` — an explicit list of what this project does and does not establish
  - `literature_review.md` — the published research this project is grounded in
- **`analyses/`** — the actual code, numbered in the order it runs
- **`results/`** — the figures and tables the report references

## Reproducing this

```bash
pip install -r requirements.txt
```

You'll need to place the following raw files in `data/raw/` yourself (see `docs/data_provenance.md` for exact accession numbers and download sources — not included in this repo due to size and data redistribution terms):

- `GSE213092_series_matrix` and `GSE213092_HRCRC.txt` (discovery cohort, from GEO)
- `TCGA-COAD.clinical.tsv.gz` and `TCGA-COAD.star_counts.tsv.gz` (external validation, from UCSC Xena)
- `hallmark_gene_sets.gmt` (MSigDB Hallmark collection, gene symbols)
- `gene_annotation_reference.csv` (Ensembl-ID-to-gene-symbol lookup table used for TCGA harmonization — source documented in `docs/data_provenance.md`)

Then run the scripts in `analyses/` in numbered order (01 through 18). One note: script `06_permutation_test_discovery_1200.py` runs all 1,200 permutations in a single process if called with no arguments (this takes roughly 15 minutes) — see the comment at the top of that script if you want to split the work into smaller chunks instead.

## The one-paragraph version, if you're skimming

I initially got a near-perfect classifier separating early- and late-onset colorectal cancer patients by gene expression. That result turned out to be a technical artifact — the patient groups had been sequenced under different lab conditions — so I rejected it and rebuilt the study on a cleaner dataset. There, a gene-level model found a modest signal that held up on an independent cohort, while a pathway-based approach and a previously published gene panel did not. The genes driving the working model point toward the tumor's surrounding tissue more than the tumor cells themselves, which is an honest and important caveat given the type of data involved. The full reasoning, including the mistake and how it was caught, is documented throughout this repository rather than tidied away.
