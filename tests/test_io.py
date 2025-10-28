from scrnax import io
import numpy as np
import pandas as pd


def test_from_raw_builds_anndata():
    X = np.array([[1, 0], [0, 2]])
    obs = pd.DataFrame({"cell_id": ["A", "B"], "batch": ["b1", "b1"]})
    var = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "gene_name": ["MT-FOO", "BAR1"],
        }
    )

    adata = io.from_raw(
        X,
        obs=obs,
        var=var,
        cell_id_col="cell_id",
        gene_id_col="gene_id",
    )

    assert adata.n_obs == 2
    assert adata.n_vars == 2
    assert "gene_name" in adata.var.columns
