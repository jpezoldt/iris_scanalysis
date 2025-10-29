from __future__ import annotations
from typing import Optional
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad


def from_mtx(
    matrix_path: str,
    cell_meta_path: str,
    gene_meta_path: str,
    cell_id_col: str = "cell_id",
    gene_id_col: str = "gene_id",
    gene_name_col: Optional[str] = "gene_name",
) -> ad.AnnData:
    """
    Create AnnData from local Matrix Market + metadata tables.

    Assumptions:
    - matrix is cells x genes.
    - cell_meta_path: TSV/CSV with one row per cell, same order as matrix rows.
    - gene_meta_path: TSV/CSV with one row per gene, same order as matrix cols.

    We do NOT mutate adata here (no log_step).
    """
    X = scipy.io.mmread(matrix_path)
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    obs = pd.read_csv(cell_meta_path, sep=None, engine="python")
    var = pd.read_csv(gene_meta_path, sep=None, engine="python")

    obs = obs.set_index(cell_id_col)
    var = var.set_index(gene_id_col)

    if gene_name_col and gene_name_col in var.columns:
        var["gene_name"] = var[gene_name_col]

    adata = ad.AnnData(X=X, obs=obs, var=var)
    return adata


def from_raw(
    X,
    obs: pd.DataFrame,
    var: pd.DataFrame,
    cell_id_col: str | None = None,
    gene_id_col: str | None = None,
) -> ad.AnnData:
    """
    Convenience for tests / power users: build AnnData from in-memory objects.

    If cell_id_col / gene_id_col are provided, use those columns as index.
    """
    obs_local = obs.copy()
    var_local = var.copy()

    if cell_id_col is not None:
        obs_local = obs_local.set_index(cell_id_col)
    if gene_id_col is not None:
        var_local = var_local.set_index(gene_id_col)

    adata = ad.AnnData(X=X, obs=obs_local, var=var_local)
    return adata
