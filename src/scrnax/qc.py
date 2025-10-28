from __future__ import annotations
import numpy as np
import anndata as ad
from .utils import log_step


def compute_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    """
    Compute QC metrics and store them in adata.obs:
      - total_counts
      - n_genes_by_counts
      - pct_counts_mito

    Also infer mitochondrial genes using adata.var['gene_name'] if present.
    """
    X = adata.X
    if hasattr(X, "toarray"):
        X_local = X.toarray()
    else:
        X_local = np.asarray(X)

    total_counts = X_local.sum(axis=1)
    n_genes_by_counts = (X_local > 0).sum(axis=1)

    if "gene_name" in adata.var.columns:
        mito_mask = adata.var["gene_name"].astype(str).str.startswith(
            ("MT-", "mt-", "Mt_")
        )
    else:
        mito_mask = adata.var.index.to_series().astype(str).str.startswith(
            ("MT-", "mt-", "Mt_")
        )

    adata.var["is_mito"] = mito_mask.values

    mito_counts = X_local[:, mito_mask.values].sum(axis=1)
    pct_counts_mito = mito_counts / total_counts * 100.0

    adata.obs["total_counts"] = total_counts
    adata.obs["n_genes_by_counts"] = n_genes_by_counts
    adata.obs["pct_counts_mito"] = pct_counts_mito

    log_step(adata, "compute_qc_metrics", {})
    return adata


def filter_cells(
    adata: ad.AnnData,
    min_counts: int = 500,
    max_counts: int = 50000,
    max_pct_mito: float = 20.0,
    min_genes: int = 200,
) -> ad.AnnData:
    """
    Return a COPY of adata with only high-quality cells.
    Also logs thresholds into filtered.uns["qc_params"] and pipeline_log.
    """
    qc_mask = (
        (adata.obs["total_counts"] >= min_counts)
        & (adata.obs["total_counts"] <= max_counts)
        & (adata.obs["pct_counts_mito"] <= max_pct_mito)
        & (adata.obs["n_genes_by_counts"] >= min_genes)
    )

    filtered = adata[qc_mask].copy()
    filtered.uns["qc_params"] = {
        "min_counts": min_counts,
        "max_counts": max_counts,
        "max_pct_mito": max_pct_mito,
        "min_genes": min_genes,
    }

    log_step(filtered, "filter_cells", filtered.uns["qc_params"])
    return filtered
