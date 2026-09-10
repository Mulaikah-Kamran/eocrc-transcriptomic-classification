# Data Provenance

## Discovery cohort: GSE213092

- **Accession:** GSE213092 (NCBI GEO)
- **BioProject:** PRJNA879084
- **Title:** "Early onset colorectal cancer tissues and late onset colorectal cancer tissues" (Ha et al., Korea Research Institute of Bioscience and Biotechnology)
- **Platform:** Illumina HiSeq 2500, single BioProject, single sequencing run context
- **Design:** EOCRC (age < 50) vs. LOCRC (age > 70) tumor tissue, colon/rectum, no matched normal
- **Initial size:** 99 samples (49 EOCRC, 50 LOCRC)
- **QC actions taken:**
  - 10 malformed rows in the deposited counts file (`GSE213092_HRCRC.txt`) dropped — these contained leftover genomic-coordinate strings instead of valid gene symbols, a data-quality artifact of the original deposit, not introduced by this analysis.
  - 1 library-size outlier removed: sample `HRCRC3109` (EOCRC), 186.8M reads vs. all other samples under 43M reads (z-score ≈ 23 relative to the rest of the cohort).
- **Final discovery cohort: 98 patients (48 EOCRC, 50 LOCRC), 33,115 genes** (gene symbols, after malformed-row removal).
- **Download:** Series matrix and supplementary counts file (`GSE213092_HRCRC.txt.gz`) via the GEO accession page (`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213092`).

## External validation cohort: TCGA-COAD

- **Source:** UCSC Xena, GDC hub (`https://xenabrowser.net/datapages/`, cohort "GDC TCGA Colon Cancer (COAD)")
- **Files used:** `TCGA-COAD.clinical.tsv.gz` (phenotype/clinical, n=562 records) and `TCGA-COAD.star_counts.tsv.gz` (gene expression, n=514 samples)
- **Data-format issue caught during audit:** the "STAR - Counts" file is actually **log2(raw_count + 1)** transformed, not raw counts as labeled — verified by reversing the transform (`2^x - 1`) and confirming exact-integer recovery on spot-checked values. Corrected before any downstream use.
- **Filtering to final cohort:**
  - 514 samples had both clinical and expression data (full overlap).
  - 471 were Primary Tumor samples; deduplicated to one sample per patient (13 patients had 2 Primary Tumor samples — resolved by keeping the most recent by `days_to_collection`), giving 458 unique patients.
  - 4 patients had missing/unusable age data and were dropped, giving the **final external validation cohort: 454 patients (53 EOCRC <50, 401 LOCRC ≥50)**.
- **Gene ID handling:** Ensembl IDs (version-stripped, e.g. `ENSG00000000003.15` → `ENSG00000000003`) mapped to gene symbols via a reference vocabulary table sourced from the BulkFormer model repository's Zenodo record (DOI 10.5281/zenodo.15744294) — used here purely as a gene-annotation lookup table, unrelated to the foundation-model investigation itself.

## Exact raw file placement for reproduction

To rerun the pipeline (`analyses/01` through `18`), place these files in `data/raw/` with these exact names:

| File | Source |
|---|---|
| `GSE213092_series_matrix` | GEO accession GSE213092, series matrix download |
| `GSE213092_HRCRC.txt` | GEO accession GSE213092, supplementary counts file |
| `TCGA-COAD.clinical.tsv.gz` | UCSC Xena, "GDC TCGA Colon Cancer (COAD)" cohort, phenotype dataset |
| `TCGA-COAD.star_counts.tsv.gz` | UCSC Xena, same cohort, "gene expression RNAseq — STAR - Counts" dataset |
| `hallmark_gene_sets.gmt` | MSigDB Human Collections, "H: hallmark gene sets", Gene Symbols format |
| `gene_annotation_reference.csv` | BulkFormer Zenodo record (DOI above), file `bulkformer_gene_info.csv` — renamed here since it's used purely as an ID-mapping table, not connected to the foundation-model investigation |

## Rejected dataset: GSE251845 + GSE196006 (the "Marx pair")

- GSE251845 (LOCRC, n=22) and GSE196006 (EOCRC, n=21) — both from Marx et al. 2024 (*Frontiers in Oncology*), a matched, published EOCRC/LOCRC study.
- **Rejected as primary discovery cohort**: cohort (age group) was completely confounded with sequencing platform/protocol — GSE251845 sequenced on NovaSeq 6000 with polyA-selected RNA, GSE196006 on MiSeq with total RNA. Every sample in each cohort shares exactly one technical setup.
- **Evidence for rejection:** a classifier trained on this pair achieved ~100% accuracy — a red flag given the literature's own finding (Marx et al.'s own PCA) that no large-scale unsupervised separation exists between EOCRC/LOCRC. Direct inspection found 6,195+ genes (≈10% of the transcriptome) completely silent in one cohort and clearly detected in the other — a signature of differing library-prep chemistry, not biology. Even after removing these genes, near-perfect separation persisted using ~13,000 "robustly expressed" genes, confirming the confound was pervasive across the transcriptome, not isolated to a subset.
- **Retained only as a secondary, explicitly labeled historical comparison** — the source of the externally nominated 8-gene panel (see below), never used as primary evidence.

## Externally nominated gene panel: Marx et al. 8 genes

ALDOB, FBXL16, IL1RN, MSLN, RAC3, SLC38A11, WBSCR27, WNT11 — nominated by Marx et al. 2024 as genes distinguishing EOCRC from LOCRC in their paired tumor-vs-normal contrast analysis. **Important terminology note:** the original paper describes these genes as the basis for a survival-predictive score, not an EOCRC-vs-LOCRC classifier. In this project, only the gene *identities* are treated as externally nominated; the classifier trained on them is newly fit here, not borrowed from the original paper.

## Excluded candidates (dataset search, documented for completeness)

| Dataset | Reason excluded |
|---|---|
| ICGC-ARGO | Scientifically attractive (975 CRC samples, 234 EOCRC) but controlled-access (DACO approval), impractical for project timeline |
| GSE175433 | NanoString (~730 immune genes), not genome-wide bulk RNA-seq |
| GSE281413 | GeoMx spatial transcriptomics, only 3 EOCRC + 3 comparator samples |
| GSE39582 | Affymetrix microarray, not RNA-seq |
| Newer Korean WGS cohort | Same underlying cohort as GSE213092, not independent |
| TCGA-READ (added to COAD) | Documented biological heterogeneity between colon and rectal tumors (differing immune infiltration, driver-mutation spectrum); serious published analyses combining them apply formal batch correction (ComBat/sva), which was judged an unnecessary additional complexity given TCGA-COAD alone already provided a statistically informative, independent validation result |
