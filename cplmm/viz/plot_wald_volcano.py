import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_wald_volcano(
    wald_df,
    pval_col: str = "P-value 1",
    fdr_col: str = "Adjusted P-value (FDR)_1",
    beta_before_col: str = "Beta 1",
    beta_after_col: str = "Beta 3",
    protein_col: str = "Protein",
    pval_threshold: float = 0.05,
    fdr_threshold: float = 0.05,
    annotate: bool = True,
    annotate_list: list = None,
    style_config: dict = None,
    export: bool = False,
    export_dir: str = "Figures",
    export_name: str = "wald_volcano",
    export_formats: list = ["pdf", "svg"]
):
    """
    Generate a volcano plot for Wald test results.

    Parameters
    ----------
    wald_df : pd.DataFrame
        DataFrame containing Wald test results.
    pval_col : str
        Column name for raw p-values.
    fdr_col : str
        Column name for FDR-adjusted p-values.
    beta_before_col : str
        Column name for slope before onset.
    beta_after_col : str
        Column name for slope after onset.
    protein_col : str
        Column name for protein names.
    pval_threshold : float
        Raw p-value significance threshold.
    fdr_threshold : float
        Adjusted p-value (FDR) significance threshold.
    annotate : bool
        Whether to annotate significant proteins.
    annotate_list : list
        Optional list of proteins to annotate (overrides automatic).
    style_config : dict
        Custom matplotlib rcParams (font size, grid, etc.).
    export : bool
        Whether to save the plot.
    export_dir : str
        Directory to save plots.
    export_name : str
        Base name for saved file.
    export_formats : list
        List of export formats (e.g., ["pdf", "svg"]).
    """
    import pandas as pd
    
    # Style
    if style_config:
        plt.rcParams.update(style_config)
    else:
        plt.rcParams.update({
            'font.family': 'serif',
            'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
            'font.size': 10,
            'axes.labelsize': 10,
            'axes.titlesize': 10,
            'legend.fontsize': 10,
            'lines.linewidth': 1.5,
            'axes.linewidth': 1.5,
            'axes.grid': True,
            'axes.grid.axis': 'y',
            'grid.color': '#DDDDDD',
            'grid.linewidth': 0.5
        })

    # Copy data to avoid modification
    df = wald_df.copy()

    # Compute metrics
    df["-log10(pval)"] = -np.log10(df[pval_col])
    df["FoldChange"] = df[beta_after_col] / df[beta_before_col]
    df["log2FoldChange"] = np.sign(df["FoldChange"]) * np.log2(np.abs(df["FoldChange"]))

    # Significance assignment
    df["Significance"] = "Not Significant"
    df.loc[(df[pval_col] < pval_threshold) & (df[fdr_col] < fdr_threshold), "Significance"] = "Significant"

    # Plot
    plt.figure(figsize=(7, 5), dpi=300)
    sns.scatterplot(
        data=df,
        x="log2FoldChange",
        y="-log10(pval)",
        hue="Significance",
        palette={"Not Significant": "gray", "Significant": "#B24745"},
        alpha=0.7,
        edgecolor="black",
        linewidth=0.5
    )

    # Threshold lines
    plt.axhline(-np.log10(pval_threshold), linestyle="dashed", color="black", linewidth=1)
    plt.axvline(-1, linestyle="dashed", color="blue", linewidth=1)
    plt.axvline(1, linestyle="dashed", color="blue", linewidth=1)

    # Labels
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(p-value)")
    plt.title("Volcano Plot for Wald Test")

    # Black spines
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(1)
        spine.set_color("black")

    plt.legend(title="Significance", loc="upper left", bbox_to_anchor=(1.02, 1), frameon=True)

    # Annotation logic
    if annotate:
        if annotate_list:
            annot_df = df[df[protein_col].isin(annotate_list)]
        else:
            annot_df = df[df["Significance"] == "Significant"]

        for _, row in annot_df.iterrows():
            plt.text(row["log2FoldChange"], row["-log10(pval)"], row[protein_col],
                     fontsize=9, ha="right", va="bottom", color="black")

    plt.tight_layout()

    # Export
    if export:
        os.makedirs(export_dir, exist_ok=True)
        for fmt in export_formats:
            plt.savefig(os.path.join(export_dir, f"{export_name}.{fmt}"), format=fmt, bbox_inches="tight")

    plt.show()
