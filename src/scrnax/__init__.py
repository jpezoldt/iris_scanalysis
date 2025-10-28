"""
scrnax: opinionated scRNA-seq preprocessing, QC, and clustering.
"""

from . import io, qc, pp, tl, vis  # re-export modules

__all__ = ["io", "qc", "pp", "tl", "vis"]

__version__ = "0.0.1"
