"""Shared utility for reading MSigDB-format .gmt gene set files.

A .gmt file has one gene set per line: name, description, then tab-separated member genes.
This is used by both the Representation B pathway pipeline and the enrichment analysis,
so it lives here once instead of being copy-pasted into both places.
"""


def parse_gmt(path, restrict_to=None):
    """Parse a .gmt file into a dict of {gene_set_name: set_of_genes}.

    If restrict_to is given (e.g. the set of genes present in our expression data),
    each gene set is intersected with it -- this keeps every downstream calculation
    working only with genes we actually have data for.
    """
    genesets = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            name, genes = parts[0], parts[2:]
            gene_set = set(genes)
            if restrict_to is not None:
                gene_set = gene_set & set(restrict_to)
            genesets[name] = gene_set
    return genesets
