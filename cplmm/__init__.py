# cplmm/__init__.py
__version__ = "0.1.0"

# -------- Preprocessing --------
from .preprocessing.preprocessing import *               # status correction, onset age
from .preprocessing.expression_prep import prepare_combined_expression

# -------- Models --------
from .models.cplmm import *                              # CPLMM fitting and slope extraction

# -------- Stats --------
from .stats.wald_test import *                           # Wald test & p-values
from .stats.mannwhitney import compare_groups_mannwhitney

# -------- Survival --------
from .survival.event_computation import compute_event_df
from .survival.km_plot import plot_km_with_threshold

# -------- Visualization --------
from .viz.plot_cplmm import plot_cplmm
from .viz.plot_wald_volcano import plot_wald_volcano
from .viz.plot_quadrant_beta import plot_quadrant_beta
from .viz.plot_pathway_bubble import plot_pathway_bubble
from .viz.plot_pathway_gene_heatmap import plot_pathway_gene_heatmap
from .viz.plot_top_pathways_bar import plot_top_pathways_bar
from .viz.plot_expression_boxplot import plot_expression_boxplot

# -------- Utils --------
from .utils.palettes import JAMA_PALETTE

__all__ = [
    # Preprocessing
    "prepare_combined_expression",

    # Models
    "fit_cplmm", "extract_slopes",  # (Assuming these are defined in cplmm.py)

    # Stats
    "wald_test", "compare_groups_mannwhitney",

    # Survival
    "compute_event_df", "plot_km_with_threshold",

    # Visualization
    "plot_cplmm", "plot_wald_volcano", "plot_quadrant_beta",
    "plot_pathway_bubble", "plot_pathway_gene_heatmap",
    "plot_top_pathways_bar", "plot_expression_boxplot",

    # Utils
    "JAMA_PALETTE"
]