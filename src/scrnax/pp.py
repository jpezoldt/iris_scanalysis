from __future__ import annotations
import scanpy as sc
import anndata as ad
from .utils import log_step


def normalize_and_log1p(adata: ad.AnnData, target_sum: float = 1e4) -> ad.AnnData:
    """
    Normalize counts per cell to `target_sum`, then log1p transform.
    Mutates `adata`.
    """
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    log_step(adata, "normalize_and_log1p", {"target_sum": target_sum})
    return adata


def find_hvgs(adata: ad.AnnData, n_top_genes: int = 2000) -> ad.AnnData:
    """
    Identify highly variable genes and annotate adata.var['highly_variable'].
    Mutates `adata`.
    """
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor="seurat_v3",
    )

    log_step(adata, "find_hvgs", {"n_top_genes": n_top_genes})
    return adata


def scale_and_pca(adata: ad.AnnData, n_comps: int = 50) -> ad.AnnData:
    """
    Scale data and run PCA.
    Stores PCA in adata.obsm['X_pca'].
    Mutates `adata`.
    """
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_comps)

    log_step(adata, "scale_and_pca", {"n_comps": n_comps})
    return adata
