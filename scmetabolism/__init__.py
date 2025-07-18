"""
scMetabolism: Quantifying the single-cell metabolism activity

A Python package for quantifying metabolism activity at the single-cell resolution.
"""

__version__ = "0.2.1"
__author__ = "Qiang Gao, Yingcheng Wu"
__email__ = "gaoqiang@fudan.edu.cn"

from .core import ScMetabolism
from .visualization import MetabolismPlotter
from .utils import load_gene_sets, preprocess_data
from .go_analysis import GOAnalysis, download_go_gene_sets

__all__ = [
    "ScMetabolism",
    "MetabolismPlotter", 
    "load_gene_sets",
    "preprocess_data",
    "GOAnalysis",
    "download_go_gene_sets"
]