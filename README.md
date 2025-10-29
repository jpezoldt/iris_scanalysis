# iscrnax

`iscrnax` is a small, structured scRNA-seq analysis toolkit built on top of `scanpy`.

The goals:
- Give researchers a **clean starting point** to run QC → filtering → normalization → HVGs → PCA → neighbors/UMAP/leiden.
- Enforce **good habits**: provenance logging, parameter capture, testable code.
- Be easy to extend in a controlled way by junior team members.

## Typical usage

```python
import iscrnax as isx

adata = isx.io.from_mtx(
    matrix_path="counts.mtx",
    cell_meta_path="cells.tsv",
    gene_meta_path="genes.tsv",
    cell_id_col="cell_id",
    gene_id_col="gene_id",
    gene_name_col="gene_name",
)

adata = isx.qc.compute_qc_metrics(adata)
adata_filt = isx.qc.filter_cells(
    adata,
    min_counts=500,
    max_counts=50000,
    max_pct_mito=20,
    min_genes=200,
)

adata_filt = isx.pp.normalize_and_log1p(adata_filt, target_sum=1e4)
adata_filt = isx.pp.find_hvgs(adata_filt, n_top_genes=2000)
adata_filt = isx.pp.scale_and_pca(adata_filt, n_comps=50)

adata_filt = isx.tl.neighbors_umap_leiden(
    adata_filt,
    n_neighbors=15,
    resolution=1.0,
)

fig = isx.vis.embedding_scatter(adata_filt, basis="umap", color="leiden")
fig.show()
```

## Provenance

Every function that mutates the AnnData object records:
- which function ran
- with which parameters

in `adata.uns["pipeline_log"]`.

## Developer rules

1. No new analysis logic lives in notebooks.
2. Every function that modifies `adata` must update `pipeline_log`.
3. Add tests in `tests/` for any new logic.
4. `io` is read-only and never overwrites user data.

## Get started
1. Clone repo: https://github.com/jpezoldt/iris_scanalysis
2. Setup virtual env in /home: `python3 -m venv .venv`
3. Start virtual env: `source /xx/.venv/bin/activate`
	Qick add bashrc: 
		`cd /home/pezoldt/software/github/iris_scanalysis`
		`activate_iris() {
source /home/pezoldt/software/github/.venvs/iris_scanalysis/bin/activate
}`
		`activate_iris`
4. Install required packages: `pip install -r requirements.txt`
5. Reflection of installation in import package: `pip install -e .[dev]`
6. One line checklist for env: # Setup (Python =3.10): create & activate venv, upgrade pip, install package in editable mode
python3.11 -m venv ~/.venvs/iris_scanalysis && source ~/.venvs/iris_scanalysis/bin/activate && python -m pip install --upgrade pip setuptools wheel && cd /home/pezoldt/software/github/iris_scanalysis && pip install -e ".[dev]"
