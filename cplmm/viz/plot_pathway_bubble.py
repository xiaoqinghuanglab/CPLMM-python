import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

def plot_pathway_bubble(
    df: pd.DataFrame,
    pathway_col: str = "Cleaned_Pathway",
    category_col: str = "BioCategory_Manual",
    source_col: str = "Source",
    logq_col: str = "LogQValue",
    gene_col: str = "Gene",
    size_scale: float = 15,
    cmap: str = "RdBu_r",
    title: str = "Pathway Enrichment by Source",
    style_config: dict = None,
    export: bool = False,
    export_dir: str = "Figures",
    export_name: str = "pathway_bubble_plot",
    export_formats: list = ["pdf", "svg"]
):
    """
    Create a compact bubble plot of pathway enrichment colored by -log10(FDR) 
    and sized by gene count.

    Parameters
    ----------
    df : pd.DataFrame
        Cleaned/categorized pathway dataframe.
    pathway_col : str
        Column name for cleaned pathway names.
    category_col : str
        Column for biological categories.
    source_col : str
        Column indicating enrichment source (DAVID, Reactome, Metascape, etc.).
    logq_col : str
        Column for log-transformed significance metric (e.g., -log10(FDR)).
    gene_col : str
        Column for gene identifiers used to count unique genes.
    size_scale : float
        Scaling factor for bubble sizes.
    cmap : str
        Colormap for significance values.
    title : str
        Plot title.
    style_config : dict
        Optional Matplotlib rcParams override.
    export : bool
        Whether to export the figure.
    export_dir : str
        Directory for figure export.
    export_name : str
        File name (without extension).
    export_formats : list
        List of formats (pdf, svg, png).
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

    # Filter valid pathways
    df_valid = df[df[pathway_col].notna()].copy()

    # ------------------ Positioning ------------------
    pathway_order_df = (
        df_valid[[pathway_col, category_col]]
        .drop_duplicates()
        .sort_values([category_col, pathway_col])
        .reset_index(drop=True)
    )

    # Assign compact positions grouped by category
    compact_gap = 1
    compact_category_offsets = {}
    compact_offset = 0

    for cat, group in pathway_order_df.groupby(category_col):
        compact_category_offsets[cat] = compact_offset
        compact_offset += len(group) + compact_gap

    pathway_order_df['compact_y_pos'] = (
        pathway_order_df.groupby(category_col).cumcount() +
        pathway_order_df[category_col].map(compact_category_offsets)
    )

    # Map y positions to main dataframe
    ypos_map = dict(zip(pathway_order_df[pathway_col], pathway_order_df['compact_y_pos']))
    df_valid['compact_y_pos'] = df_valid[pathway_col].map(ypos_map)

    # Collapse gene info to count per pathway/source
    plot_df = (
        df_valid.groupby([pathway_col, source_col, 'compact_y_pos'])
        .agg(LogQValue=(logq_col, 'max'),
             GeneCount=(gene_col, 'nunique'))
        .reset_index()
    )

    # Map x positions for sources
    sources = plot_df[source_col].unique()
    x_map = dict(zip(sources, range(len(sources))))
    plot_df['x_pos'] = plot_df[source_col].map(x_map)

    # ------------------ Plot ------------------
    fig, ax = plt.subplots(figsize=(7.5, (compact_offset + 2) * 0.5), dpi=300)

    scatter = ax.scatter(
        x=plot_df['x_pos'],
        y=plot_df['compact_y_pos'],
        s=plot_df['GeneCount'] * size_scale,
        c=plot_df['LogQValue'],
        cmap=cmap,
        edgecolor='black'
    )

    # Y-axis: pathways
    ax.set_yticks(pathway_order_df['compact_y_pos'])
    ax.set_yticklabels(pathway_order_df[pathway_col])
    ax.invert_yaxis()

    # X-axis: sources
    ax.set_xticks(range(len(sources)))
    ax.set_xticklabels(sources, rotation=45, ha='right')

    # Category separators
    category_ends = pathway_order_df.groupby(category_col)['compact_y_pos'].max()
    for end in category_ends:
        ax.axhline(y=end + 0.5, color='gray', linestyle='--', linewidth=0.6)

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.5)
    cbar.set_label('-log10(FDR)', rotation=270, labelpad=12)
    cbar.ax.tick_params(labelsize=10)

    # Bubble size legend (quartile sizes)
    unique_sizes = plot_df['GeneCount'].dropna().unique()
    legend_sizes = np.percentile(unique_sizes, [25, 50, 75])
    legend_sizes = np.unique(np.round(legend_sizes).astype(int))

    for size in legend_sizes:
        ax.scatter([], [], s=size * size_scale, c='gray', edgecolors='black', label=f'{size} genes')

    legend = ax.legend(
        title="Gene Count",
        loc='center left',
        bbox_to_anchor=(1.05, 0.7),
        frameon=True,
        edgecolor='black'
    )
    legend.get_title().set_fontsize(12)
    legend.get_frame().set_linewidth(1.0)

    # Spines & title
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color('black')

    plt.title(title, fontsize=16, pad=10)
    plt.tight_layout()

    # Export
    if export:
        os.makedirs(export_dir, exist_ok=True)
        for fmt in export_formats:
            plt.savefig(os.path.join(export_dir, f"{export_name}.{fmt}"),
                        format=fmt, bbox_inches="tight")

    plt.show()
