from scrnax import pp, qc


def test_normalize_and_log1p_runs(tiny_adata):
    adata = qc.compute_qc_metrics(tiny_adata.copy())
    adata = pp.normalize_and_log1p(adata)
    assert "pipeline_log" in adata.uns


def test_find_hvgs_runs(tiny_adata):
    adata = qc.compute_qc_metrics(tiny_adata.copy())
    adata = pp.normalize_and_log1p(adata)
    adata = pp.find_hvgs(adata, n_top_genes=2)
    assert "highly_variable" in adata.var.columns


def test_scale_and_pca_runs(tiny_adata):
    adata = qc.compute_qc_metrics(tiny_adata.copy())
    adata = pp.normalize_and_log1p(adata)
    adata = pp.find_hvgs(adata, n_top_genes=2)
    adata = pp.scale_and_pca(adata, n_comps=3)
    assert "X_pca" in adata.obsm
