import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
import pandas as pd
from matplotlib.lines import Line2D

def plot_pathway_gene_heatmap(
    df: pd.DataFrame,
    pathway_col: str = "Cleaned_Pathway",
    category_col: str = "BioCategory_Manual",
    gene_col: str = "Gene",
    palette: list = None,
    style_config: dict = None,
    title: str = "Pathway-Gene Membership Heatmap",
    export: bool = False,
    export_dir: str = "Figures",
    export_name: str = "pathway_gene_heatmap",
    export_formats: list = ["pdf", "svg"]
):
    """
    Create a heatmap showing gene membership across pathways, grouped by biological category.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing pathways, categories, and genes.
    pathway_col : str
        Column name for pathways.
    category_col : str
        Column name for pathway categories.
    gene_col : str
        Column name for genes.
    palette : list
        List of colors for categories (if None, a default palette is used).
    style_config : dict
        Optional matplotlib rcParams override.
    title : str
        Plot title.
    export : bool
        Whether to save the figure.
    export_dir : str
        Directory to save plots.
    export_name : str
        Filename (without extension).
    export_formats : list
        Formats to export (e.g., pdf, svg, png).
    """

    # ------------------ Style ------------------
    if style_config:
        plt.rcParams.update(style_config)
    else:
        plt.rcParams.update({
            'font.family': 'serif',
            'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
            'font.size': 14,
            'axes.labelsize': 14,
            'axes.titlesize': 16,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 12,
            'axes.linewidth': 1.5,
            'pdf.fonttype': 42
        })

    # ------------------ Prepare Data ------------------
    df_valid = df[df[pathway_col].notna()].copy()

    # Sort pathways by category and name
    pathway_order_df = (
        df_valid[[pathway_col, category_col]]
        .drop_duplicates()
        .sort_values([category_col, pathway_col])
        .reset_index(drop=True)
    )

    # Assign compact Y positions grouped by category
    gap = 1
    category_offsets = {}
    offset = 0
    for cat, group in pathway_order_df.groupby(category_col):
        category_offsets[cat] = offset
        offset += len(group) + gap

    pathway_order_df["compact_y_pos"] = (
        pathway_order_df.groupby(category_col).cumcount()
        + pathway_order_df[category_col].map(category_offsets)
    )

    # ------------------ Build Pathway-Gene Matrix ------------------
    heatmap_df = df_valid[[pathway_col, gene_col]].drop_duplicates()
    heatmap_df["Present"] = 1
    matrix = heatmap_df.pivot_table(index=pathway_col, columns=gene_col, values="Present", fill_value=0)
    matrix = matrix.reindex(pathway_order_df[pathway_col])

    # Melt to long format
    heatmap_long = matrix.reset_index().melt(id_vars=pathway_col, var_name=gene_col, value_name="Present")
    heatmap_long = heatmap_long[heatmap_long["Present"] == 1]
    heatmap_long = heatmap_long.merge(pathway_order_df, on=pathway_col, how="left")

    # Assign X positions for genes
    gene_order = sorted(heatmap_long[gene_col].unique())
    gene_to_x = {gene: i for i, gene in enumerate(gene_order)}
    heatmap_long["x_pos"] = heatmap_long[gene_col].map(gene_to_x)

    # ------------------ Palette ------------------
    bio_categories = sorted(df_valid[category_col].dropna().unique())
    if palette is None:
        default_palette = [
            '#DF8F44', '#00A1D5', '#B24745',
            '#79AF97', '#6A6599', '#374E55', '#80796B',
            '#AA4499', '#117733', '#999933', '#882255'
        ]
        palette = default_palette[:len(bio_categories)]
    category_palette = dict(zip(bio_categories, palette))

    # ------------------ Plot ------------------
    fig, ax = plt.subplots(figsize=(len(gene_order) * 0.6, (offset + 2) * 0.5), dpi=300)

    # Draw colored squares
    for _, row in heatmap_long.iterrows():
        rect = patches.Rectangle(
            (row["x_pos"] - 0.5, row["compact_y_pos"] - 0.5),
            1, 1,
            linewidth=0.3,
            edgecolor="black",
            facecolor=category_palette[row[category_col]]
        )
        ax.add_patch(rect)

    # Axes setup
    ax.set_xlim(-0.5, len(gene_order) - 0.5)
    ax.set_ylim(offset + 1, -1)
    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=90)
    ax.set_yticks(pathway_order_df["compact_y_pos"])
    ax.set_yticklabels(pathway_order_df[pathway_col])

    # Grid lines
    for x in range(len(gene_order) + 1):
        ax.axvline(x - 0.5, color="gray", linewidth=0.3, zorder=1)
    for y in range(offset + 2):
        ax.axhline(y - 0.5, color="gray", linewidth=0.3, zorder=1)

    # Category dividers
    for end in pathway_order_df.groupby(category_col)["compact_y_pos"].max():
        ax.axhline(y=end + 0.5, color="gray", linestyle="--", linewidth=0.6)

    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color("black")

    # Legend
    legend_elements = [
        Line2D([0], [0], marker="s", color="w", label=cat,
               markerfacecolor=category_palette[cat], markeredgecolor="black", markersize=10)
        for cat in bio_categories
    ]
    legend = ax.legend(
        handles=legend_elements,
        title="Category",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=True
    )
    legend.get_title().set_fontsize(12)
    legend.get_frame().set_linewidth(1.0)

    # Title
    plt.title(title, fontsize=16, pad=10)

    # Export
    plt.tight_layout()
    if export:
        os.makedirs(export_dir, exist_ok=True)
        for fmt in export_formats:
            plt.savefig(os.path.join(export_dir, f"{export_name}.{fmt}"), format=fmt, bbox_inches="tight")

    plt.show()
