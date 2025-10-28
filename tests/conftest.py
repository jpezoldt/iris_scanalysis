import pytest
import numpy as np
import pandas as pd
import anndata as ad


@pytest.fixture
def tiny_adata():
    """
    Minimal AnnData for tests:
    3 cells x 4 genes, with mitochondrial-looking gene names.
    """
    X = np.array([
        [10, 0, 3, 1],
        [5,  2, 0, 0],
        [0,  1, 0, 7],
    ])

    obs = pd.DataFrame(
        {
            "cell_id": ["cellA", "cellB", "cellC"],
        }
    ).set_index("cell_id")

    var = pd.DataFrame(
        {
            "gene_id":   ["g1",     "g2",    "g3",    "g4"],
            "gene_name": ["MT-CO1", "GENE2", "GENE3", "MT-ND1"],
        }
    ).set_index("gene_id")

    adata = ad.AnnData(X=X, obs=obs, var=var)
    return adata
