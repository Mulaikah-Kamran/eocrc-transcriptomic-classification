"""Shared utilities for preparing TCGA-COAD data to match the GSE213092 gene/sample space.

Centralizing this here (rather than duplicating it in every script that touches TCGA)
is a direct fix for a real problem: two earlier versions of this pipeline re-derived
the same TCGA processing from scratch instead of reusing one shared, tested version.
"""
import pandas as pd
import numpy as np


def reverse_xena_log2_transform(log2_df):
    """UCSC Xena's 'STAR - Counts' file for TCGA-COAD is actually log2(raw_count + 1),
    not raw counts as labeled -- caught during the original data audit by reversing
    the transform and confirming it recovers clean integer read counts. This function
    does that reversal so downstream code can treat the result as raw counts.
    """
    return (2 ** log2_df) - 1


def map_ensembl_to_symbol(raw_counts_df, gene_info_df):
    """Strip Ensembl version suffixes (e.g. '.15') and collapse to gene symbols
    using a reference annotation table (ensg_id <-> gene_symbol columns expected).
    Multiple Ensembl IDs mapping to the same symbol are summed.
    """
    df = raw_counts_df.copy()
    df.index = df.index.str.split(".").str[0]
    ensg_to_symbol = (
        gene_info_df.drop_duplicates(subset="ensg_id", keep=False)
        .set_index("ensg_id")["gene_symbol"]
        .to_dict()
    )
    df["symbol"] = df.index.map(lambda g: ensg_to_symbol.get(g))
    df = df.dropna(subset=["symbol"]).groupby("symbol").sum()
    return df


def counts_to_log2cpm(counts_df):
    """Standard log2(CPM + 1) normalization, used identically for GSE213092 and TCGA
    so both cohorts are put into a genuinely comparable representation space."""
    libsize = counts_df.sum(axis=0)
    cpm = counts_df.div(libsize, axis=1) * 1e6
    return np.log2(cpm + 1)


def apply_verified_renames(log2cpm_df, renames_dict):
    """Apply only manually-verified gene-symbol renames (e.g. IL8 -> CXCL8), so a
    frozen model's original gene name can still find its data under a newer name.
    See docs/decision_log.md for why only verified renames are applied here, not
    a blanket automated renaming pass.
    """
    df = log2cpm_df.copy()
    for old_name, new_name in renames_dict.items():
        if new_name in df.index and old_name not in df.index:
            df.loc[old_name] = df.loc[new_name]
    return df
