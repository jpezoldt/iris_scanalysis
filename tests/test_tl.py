from scrnax import pp, qc, tl


def test_neighbors_umap_leiden_runs(tiny_adata):
    adata = qc.compute_qc_metrics(tiny_adata.copy())
    adata = pp.normalize_and_log1p(adata)
    adata = pp.find_hvgs(adata, n_top_genes=2)
    adata = pp.scale_and_pca(adata, n_comps=3)

    adata = tl.neighbors_umap_leiden(
        adata,
        n_neighbors=2,
        resolution=0.5,
    )

    assert "X_umap" in adata.obsm
    assert "leiden" in adata.obs
    assert "pipeline_log" in adata.uns
