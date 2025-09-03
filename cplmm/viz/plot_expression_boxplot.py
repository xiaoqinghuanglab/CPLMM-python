import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_expression_boxplot(
    combined_expr,
    gene_order,
    hue_col="Source",
    expression_col="Expression",
    gene_col="Gene",
    palette=None,
    title="Expression of Top Genes Across Categories",
    y_limit=None,
    figsize=(32, 8),
    rotation=45,
    export=False,
    export_dir="Figures",
    export_name="gene_expression_boxplot",
    export_formats=["pdf", "svg"]
):
    """
    Plot expression of genes across diagnostic categories as boxplots.

    Parameters
    ----------
    combined_expr : pd.DataFrame
        Long-format dataframe with columns: Gene, Expression, Source.
    gene_order : list
        Ordered list of genes to plot on the x-axis.
    hue_col : str
        Column for grouping (e.g., diagnostic categories).
    expression_col : str
        Column containing expression values.
    gene_col : str
        Column containing gene names.
    palette : dict
        Color palette mapping categories to colors.
    title : str
        Plot title.
    y_limit : float or None
        Upper limit for y-axis (default None = auto-scale).
    figsize : tuple
        Figure size.
    rotation : int
        X-axis label rotation.
    export : bool
        Whether to export the plot.
    export_dir : str
        Directory for saving plots.
    export_name : str
        Base filename for exported plots.
    export_formats : list
        Formats for export (e.g., pdf, svg, png).
    """
    # Ensure gene order is categorical
    combined_expr[gene_col] = pd.Categorical(combined_expr[gene_col], categories=gene_order, ordered=True)

    # Default JAMA palette if not provided
    if palette is None:
        palette = {
            "Normal_only": "#00A1D5",
            "SCD": "#DF8F44",
            "MCI": "#B24745",
            "AD_Dementia": "#79AF97",
            "FTD_Dementia": "#6A6599"
        }

    # Set style
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size": 20,
        "axes.labelsize": 20,
        "axes.titlesize": 22,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 18,
        "axes.linewidth": 1.5,
        "pdf.fonttype": 42
    })

    # Flier properties (outliers)
    flier_props = dict(
        marker="o", markerfacecolor="white",
        markeredgecolor="black", markersize=2,
        markeredgewidth=0.5, linestyle="none"
    )

    # Plot
    fig, ax = plt.subplots(figsize=figsize, dpi=300)
    sns.boxplot(
        data=combined_expr,
        x=gene_col,
        y=expression_col,
        hue=hue_col,
        palette=palette,
        linewidth=0.8,
        width=0.9,
        flierprops=flier_props,
        ax=ax
    )

    # Aesthetics
    if y_limit:
        ax.set_ylim(top=y_limit)
    ax.set_title(title, pad=10)
    ax.set_xlabel("Gene")
    ax.set_ylabel("Expression")
    plt.xticks(rotation=rotation, ha="right")

    # Legend
    legend = ax.legend(title="Category", bbox_to_anchor=(1.01, 1), loc="upper left", frameon=True)
    legend.get_frame().set_linewidth(1.0)
    legend.get_frame().set_edgecolor("black")

    # Grid and spines
    ax.yaxis.grid(True, linewidth=0.4)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color("black")

    plt.tight_layout()

    # Export
    if export:
        os.makedirs(export_dir, exist_ok=True)
        for fmt in export_formats:
            plt.savefig(f"{export_dir}/{export_name}.{fmt}", format=fmt, bbox_inches="tight")

    plt.show()
