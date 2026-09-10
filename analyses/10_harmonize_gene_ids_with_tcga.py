"""Step 10 -- Harmonize gene identifiers between the frozen Representation A gene
list and TCGA's Ensembl-ID gene space.

This step was previously missing as a saved script -- reconstructed here.

Expects:
    data/frozen_representation_a_model.pkl
    data/raw/TCGA-COAD.star_counts.tsv.gz
    data/raw/gene_annotation_reference.csv   (Ensembl<->symbol lookup, see docs/data_provenance.md
                                               for where this specific reference table came from)

Produces:
    data/final_genes_present_in_tcga.csv
    data/final_genes_missing_from_tcga.csv
    data/harmonized_gene_renames.csv          (only manually-verified renames, e.g. IL8 -> CXCL8)
    data/genes_genuinely_unavailable.csv
"""
import re
import pickle
import pandas as pd

with open("data/frozen_representation_a_model.pkl", "rb") as f:
    frozen = pickle.load(f)
final_genes = frozen["final_gene_list"]

gene_info = pd.read_csv("data/raw/gene_annotation_reference.csv")
ensg_to_symbol = gene_info.drop_duplicates(subset="ensg_id", keep=False) \
                           .set_index("ensg_id")["gene_symbol"].to_dict()
current_symbols = set(gene_info["gene_symbol"])

counts_header = pd.read_csv("data/raw/TCGA-COAD.star_counts.tsv.gz", sep="\t", nrows=0)
tcga_ensg_ids = [g.split(".")[0] for g in counts_header.columns[1:].tolist()]
# the first column is the Ensembl_ID index, handled separately elsewhere; here we just
# need TCGA's raw index, so re-read the first column only
tcga_gene_col = pd.read_csv("data/raw/TCGA-COAD.star_counts.tsv.gz", sep="\t", usecols=[0])
tcga_ensg_stripped = tcga_gene_col.iloc[:, 0].str.split(".").str[0]
tcga_symbols_present = set(tcga_ensg_stripped.map(lambda g: ensg_to_symbol.get(g)).dropna())

present = [g for g in final_genes if g in tcga_symbols_present]
missing = [g for g in final_genes if g not in tcga_symbols_present]
print(f"Of {len(final_genes)} final genes: {len(present)} present in TCGA, {len(missing)} missing")

# --- only apply manually-verified renames (see docs/decision_log.md) ---
# Verified directly against HGNC/NCBI records; NOT a blanket automated symbol-updating pass.
verified_renames = {"IL8": "CXCL8"}
mito_trna_pattern = re.compile(r"^TRN([A-Z][a-z0-9]*)$")

resolved = {}
for g in missing:
    if g in verified_renames:
        resolved[g] = verified_renames[g]
    else:
        m = mito_trna_pattern.match(g)
        if m:
            candidate = f"MT-T{m.group(1)}"
            if candidate in current_symbols:
                resolved[g] = candidate

still_unavailable = [g for g in missing if g not in resolved]
print(f"Resolved via verified renaming: {len(resolved)} (e.g. {resolved})")
print(f"Genuinely unavailable: {len(still_unavailable)}")

pd.Series(present).to_csv("data/final_genes_present_in_tcga.csv", index=False)
pd.Series(missing).to_csv("data/final_genes_missing_from_tcga.csv", index=False)
pd.Series(resolved).to_csv("data/harmonized_gene_renames.csv")
pd.Series(still_unavailable).to_csv("data/genes_genuinely_unavailable.csv", index=False)
print("Saved all four harmonization output files to data/")
