import os
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from .event_computation import compute_event_df

def plot_km_with_threshold(
    biomarker_name,
    threshold,
    wd_df,
    normal_df,
    abnm_df,
    threshold_dict=None,
    save=False,
    save_path="Figures",
    jama_palette=None,
    time_points=np.arange(-10, 6, 2)
):
    """
    Plot Kaplan-Meier curves for biomarker threshold crossing, including log-rank test and at-risk table.

    Parameters
    ----------
    biomarker_name : str
        Name of biomarker column (e.g., 'NFL').
    threshold : float
        Predefined biomarker threshold.
    wd_df : pd.DataFrame
        Status-change patients dataframe.
    normal_df : pd.DataFrame
        Always-normal patients dataframe.
    abnm_df : pd.DataFrame
        Always-abnormal patients dataframe.
    threshold_dict : dict
        Optional threshold lookup dictionary.
    save : bool
        Whether to save figure.
    save_path : str
        Directory for saving figure.
    jama_palette : dict
        Optional custom color palette for groups.
    time_points : array-like
        Time points for plotting and at-risk table.
    """
    # Extract MCI subgroup from status-change patients
    mci_subids = wd_df.loc[wd_df["CATEGORY"] == "MCI", "SUBID"].unique()
    mci_df = wd_df[wd_df["SUBID"].isin(mci_subids)]

    # Compute event data
    event_df_normal = compute_event_df(normal_df, threshold, biomarker_name, "Normal")
    event_df_mci = compute_event_df(mci_df, threshold, biomarker_name, "MCI")
    event_df_statuschange = compute_event_df(wd_df, threshold, biomarker_name, "Status-Change")
    event_df = pd.concat([event_df_statuschange, event_df_normal, event_df_mci], ignore_index=True)

    # Default palette
    if jama_palette is None:
        jama_palette = {"Status-Change": "#DF8F44", "Normal": "#00A1D5", "MCI": "#B24745"}

    # Fit KM curves
    fig, ax = plt.subplots(figsize=(6, 8), dpi=300)
    kmf_dict, at_risk_dict = {}, {}
    for label, group_df in event_df.groupby("group"):
        kmf = KaplanMeierFitter()
        kmf.fit(group_df["time"], group_df["event"], label=label)
        kmf.plot_survival_function(ax=ax, ci_show=True, linewidth=1.5, color=jama_palette.get(label, 'black'))
        kmf_dict[label] = kmf
        at_risk_dict[label] = [
            kmf.event_table.loc[kmf.event_table.index >= t, 'at_risk'].iloc[0] if any(kmf.event_table.index >= t) else 0
            for t in time_points
        ]

    # Plot labels and axis
    ax.axvline(0, linestyle='--', color='gray', linewidth=1)
    ax.set_title(f"{biomarker_name} ≥ {threshold:.2f}", fontsize=20)
    ax.set_xlabel("Years to Onset", fontsize=20)
    ax.set_ylabel("Proportion Not Yet Crossed", fontsize=20)
    ax.set_xticks(time_points)

    # Log-rank test (Normal vs MCI)
    group1 = event_df[event_df["group"] == "Normal"]
    group2 = event_df[event_df["group"] == "MCI"]
    result = logrank_test(group1["time"], group1["event"], group2["time"], group2["event"])
    ax.annotate(f"p = {result.p_value:.3e}", xy=(0.65, 0.1), xycoords='axes fraction', fontsize=18)

    # At-risk table
    plt.subplots_adjust(bottom=0.28)
    table_y, line_spacing = -0.25, 0.06
    for i, t in enumerate(time_points):
        ax.text(t, table_y, f"{t}", ha='center', va='center', fontsize=18)
    for j, (label, at_risk) in enumerate(at_risk_dict.items()):
        for i, n in enumerate(at_risk):
            ax.text(time_points[i], table_y - (j+1)*line_spacing, str(n), ha='center', va='center', fontsize=18)
        ax.text(time_points[0]-1.5, table_y - (j+1)*line_spacing, label, ha='right', va='center', fontsize=18)

    # Style and save
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color("black")
    ax.yaxis.grid(True, linestyle=':', linewidth=0.4)
    ax.set_axisbelow(True)

    if save:
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(os.path.join(save_path, f"KM_{biomarker_name}_JAMA.svg"), format="svg", bbox_inches="tight")

    plt.show()
