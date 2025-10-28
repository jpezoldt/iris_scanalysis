from __future__ import annotations
from typing import Any, Dict
import anndata as ad


def log_step(adata: ad.AnnData, step: str, params: Dict[str, Any]) -> None:
    """
    Append a provenance record to adata.uns['pipeline_log'].
    """
    if "pipeline_log" not in adata.uns:
        adata.uns["pipeline_log"] = []
    adata.uns["pipeline_log"].append({"step": step, "params": params})
