from __future__ import annotations
import scanpy as sc
import anndata as ad
from .utils import log_step


def neighbors_umap_leiden(
    adata: ad.AnnData,
    n_neighbors: int = 15,
    resolution: float = 1.0,
) -> ad.AnnData:
    """
    Build neighborhood graph, compute UMAP, run Leiden clustering.
    Adds:
      - adata.obsm['X_umap']
      - adata.obs['leiden']
    Mutates `adata`.
    """
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)

    log_step(
        adata,
        "neighbors_umap_leiden",
        {"n_neighbors": n_neighbors, "resolution": resolution},
    )
    return adata
