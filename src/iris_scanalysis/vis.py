from __future__ import annotations
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad


def embedding_scatter(
    adata: ad.AnnData,
    basis: str = "umap",
    color: str = "leiden",
):
    """
    Return a matplotlib Figure showing low-D embedding colored by `color`.
    We don't plt.show() here; caller decides.
    """
    sc.pl.embedding(
        adata,
        basis=basis,
        color=color,
        show=False,
    )
    fig = plt.gcf()
    return fig


def violin_qc(
    adata: ad.AnnData,
    keys=None,
):
    """
    Return violin plot Figure for QC metrics.
    """
    if keys is None:
        keys = ["total_counts", "n_genes_by_counts", "pct_counts_mito"]

    sc.pl.violin(
        adata,
        keys=keys,
        groupby=None,
        rotation=45,
        show=False,
    )
    fig = plt.gcf()
    return fig
