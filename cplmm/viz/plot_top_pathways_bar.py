import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def plot_top_pathways_bar(
    df: pd.DataFrame,
    pathway_col: str = "Cleaned_Pathway",
    gene_col: str = "Gene",
    logq_col: str = "LogQValue",
    category_col: str = "BioCategory_Manual",
    top_n: int = 15,
    palette: list = None,
    annotate: bool = True,
    annotation_format: str = "LogQ = {val:.2f}",
    title: str = "Top Pathways by Gene Count",
    export: bool = False,
    export_dir: str = "Figures",
    export_name: str = "top_pathways_bar",
    export_formats: list = ["pdf", "svg"]
):
    """
    Plot top N pathways ranked by gene count, annotated with log(q-values).

    Parameters
    ----------
    df : pd.DataFrame
        Pathway enrichment dataframe with pathway, gene, category, and log(q-value).
    pathway_col : str
        Column for pathway names.
    gene_col : str
        Column for gene identifiers.
    logq_col : str
        Column for log-transformed q-values.
    category_col : str
        Column for biological categories.
    top_n : int
        Number of pathways to display (default 15).
    palette : list
        Custom color palette for categories.
    annotate : bool
        Whether to annotate bars with log(q) or q-values.
    annotation_format : str
        Format for annotation (e.g., "LogQ = {val:.2f}" or "q = {val:.1e}").
    title : str
        Title of the plot.
    export : bool
        Whether to export the figure.
    export_dir : str
        Directory to save figures.
    export_name : str
        Base file name for exports.
    export_formats : list
        List of formats to export (e.g., pdf, svg, png).
    """
    
    # Prepare top pathways
    top_pathways = (
        df.groupby(pathway_col)[gene_col]
        .nunique()
        .nlargest(top_n)
        .reset_index()
        .rename(columns={gene_col: "GeneCount"})
    )

    # Merge with logQ and categories
    logqvals = df.groupby(pathway_col)[logq_col].max().reset_index()
    cats = df.groupby(pathway_col)[category_col].first().reset_index()

    top_pathways = top_pathways.merge(logqvals, on=pathway_col)
    top_pathways = top_pathways.merge(cats, on=pathway_col)
    top_pathways = top_pathways.sort_values("GeneCount", ascending=False)

    # Prepare palette
    bio_categories = sorted(df[category_col].dropna().unique())
    if palette is None:
        default_palette = [
            '#DF8F44', '#00A1D5', '#B24745',
            '#79AF97', '#6A6599', '#374E55', '#80796B',
            '#AA4499', '#117733', '#999933', '#882255'
        ]
        palette = default_palette[:len(bio_categories)]
    category_palette = dict(zip(bio_categories, palette))

    # Plot
    plt.figure(figsize=(9, 6))
    ax = sns.barplot(
        data=top_pathways,
        x="GeneCount",
        y=pathway_col,
        hue=category_col,
        dodge=False,
        palette=category_palette
    )

    # Annotate bars with log(q) values
    if annotate:
        for i, row in top_pathways.iterrows():
            ax.text(
                row["GeneCount"] + 0.3,
                list(top_pathways.index).index(i),
                annotation_format.format(val=row[logq_col]),
                va="center",
                ha="left",
                fontsize=9
            )

    # Formatting
    ax.set_xlabel("Gene Count")
    ax.set_ylabel("Pathway")
    ax.set_title(title)
    ax.legend(title="Biological Category", bbox_to_anchor=(1.02, 0.5), loc="upper left")
    plt.tight_layout()

    # Export
    if export:
        import os
        os.makedirs(export_dir, exist_ok=True)
        for fmt in export_formats:
            plt.savefig(f"{export_dir}/{export_name}.{fmt}", format=fmt, bbox_inches="tight")

    plt.show()
