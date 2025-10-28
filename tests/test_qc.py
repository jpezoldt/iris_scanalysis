from scrnax import qc


def test_compute_qc_metrics_adds_columns(tiny_adata):
    adata = qc.compute_qc_metrics(tiny_adata)
    assert "total_counts" in adata.obs.columns
    assert "n_genes_by_counts" in adata.obs.columns
    assert "pct_counts_mito" in adata.obs.columns
    assert "is_mito" in adata.var.columns
    assert "pipeline_log" in adata.uns


def test_filter_cells_returns_subset(tiny_adata):
    adata = qc.compute_qc_metrics(tiny_adata)
    filtered = qc.filter_cells(
        adata,
        min_counts=1,
        max_counts=100,
        max_pct_mito=100,
        min_genes=1,
    )

    assert filtered.n_obs <= adata.n_obs
    assert "qc_params" in filtered.uns
    assert "pipeline_log" in filtered.uns
