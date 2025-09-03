import numpy as np
import matplotlib.pyplot as plt
import os

def plot_quadrant_beta(
    wald_df,
    beta_x_col: str = "Beta 1",
    beta_y_col: str = "Beta 3",
    fdr_col: str = "Adjusted P-value (FDR)_1",
    protein_col: str = "Protein",
    fdr_threshold: float = 0.05,
    annotate: bool = True,
    style_config: dict = None,
    export: bool = False,
    export_dir: str = "Figures",
    export_name: str = "quadrant_plot",
    export_formats: list = ["pdf", "svg"]
):
    """
    Generate a quadrant plot of beta coefficients from Wald test results.

    Parameters
    ----------
    wald_df : pd.DataFrame
        DataFrame containing Wald test results.
    beta_x_col : str
        Column for x-axis beta coefficient.
    beta_y_col : str
        Column for y-axis beta coefficient.
    fdr_col : str
        Column for FDR-adjusted p-values.
    protein_col : str
        Column for protein identifiers.
    fdr_threshold : float
        Significance threshold for FDR.
    annotate : bool
        Whether to annotate significant proteins.
    style_config : dict
        Optional matplotlib rcParams.
    export : bool
        Whether to save the plot.
    export_dir : str
        Directory to save plots.
    export_name : str
        Base filename for exported plots.
    export_formats : list
        Formats to export (e.g., ["pdf", "svg", "png"]).
    """
    
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

    # Extract data
    x = wald_df[beta_x_col].values
    y = wald_df[beta_y_col].values
    fdr = wald_df[fdr_col].values
    proteins = wald_df[protein_col].values

    # Identify significant points
    significant = fdr < fdr_threshold

    # Plot
    plt.figure(figsize=(6, 6), dpi=300)

    # Non-significant points
    plt.scatter(x[~significant], y[~significant], color='black', s=10, alpha=0.3, label='Non-significant')

    # Significant points
    plt.scatter(x[significant], y[significant], color='#B24745', s=20, alpha=0.8, label=f'Significant (FDR < {fdr_threshold})')

    # Axes at zero
    plt.axhline(0, color='black', linewidth=1)
    plt.axvline(0, color='black', linewidth=1)

    # Annotate significant proteins
    if annotate:
        for i in np.where(significant)[0]:
            plt.text(x[i], y[i], proteins[i], fontsize=9, color='#B24745', ha='right', va='bottom')

    # Quadrant counts
    q1 = np.sum((x > 0) & (y > 0))
    q2 = np.sum((x < 0) & (y > 0))
    q3 = np.sum((x < 0) & (y < 0))
    q4 = np.sum((x > 0) & (y < 0))

    plt.text(max(x) * 0.7, max(y) * 0.9, f'n={q1}', fontsize=10)
    plt.text(min(x) * 0.7, max(y) * 0.9, f'n={q2}', fontsize=10)
    plt.text(min(x) * 0.7, min(y) * 0.9, f'n={q3}', fontsize=10)
    plt.text(max(x) * 0.7, min(y) * 0.9, f'n={q4}', fontsize=10)

    # Labels
    plt.xlabel(f"{beta_x_col} (Before Onset)" if "1" in beta_x_col else beta_x_col)
    plt.ylabel(f"{beta_y_col} (After Onset)" if "3" in beta_y_col else beta_y_col)
    plt.title("Quadrant Plot of Beta Coefficients")

    # Spines
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color('black')

    # Axis limits with 10% padding
    plt.xlim([min(x) * 1.1, max(x) * 1.1])
    plt.ylim([min(y) * 1.1, max(y) * 1.1])

    plt.legend(frameon=False, loc='best')
    plt.tight_layout()

    # Export
    if export:
        os.makedirs(export_dir, exist_ok=True)
        for fmt in export_formats:
            plt.savefig(os.path.join(export_dir, f"{export_name}.{fmt}"), format=fmt, bbox_inches="tight")

    plt.show()
