import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

def plot_cplmm(
    df_status_change: pd.DataFrame,
    df_normal: pd.DataFrame,
    df_abnormal: pd.DataFrame,
    protein: str,
    covariates: list[str] = ["SEX", "BASELINE_AGE"],
    years_since_onset_col: str = "years_since_onset",
    subject_id_col: str = "SUBID",
    group_colors: dict = None,
    group_labels: dict = None,
    style_config: dict = None,
    export: bool = False,
    export_dir: str = "Figures",
    export_formats: list = ["pdf", "svg"]
):
    """
    Plot change-point LMM trajectories for protein expression across patient groups.

    Parameters
    ----------
    df_status_change : pd.DataFrame
        Patients with a status change (Normal → Abnormal).
    df_normal : pd.DataFrame
        Patients who stayed Normal throughout.
    df_abnormal : pd.DataFrame
        Patients who stayed Abnormal throughout.
    protein : str
        Protein column to analyze.
    covariates : list
        Covariates to include in the model.
    years_since_onset_col : str
        Column indicating years since onset.
    subject_id_col : str
        Column identifying patient ID.
    group_colors : dict
        Colors for groups. Keys: 'status_change', 'normal', 'abnormal'.
    group_labels : dict
        Labels for groups for plotting. Same keys as above.
    style_config : dict
        Optional Matplotlib rcParams override.
    export : bool
        Whether to save plots to file.
    export_dir : str
        Directory to save figures to.
    export_formats : list
        List of formats to export (e.g. ["pdf", "svg"]).
    """

    if style_config:
        plt.rcParams.update(style_config)

    # Defaults
    group_colors = group_colors or {
        'status_change': '#DF8F44',
        'normal': '#00A1D5',
        'abnormal': '#B24745'
    }

    group_labels = group_labels or {
        'status_change': 'Status Change',
        'normal': 'Normal',
        'abnormal': 'Abnormal'
    }

    mean_covariates = {
        cov: df_status_change[cov].mean() for cov in covariates if cov in df_status_change.columns
    }

    def fit_and_plot(df, formula_rhs, group_key):
        if df.empty:
            return

        df = df.copy()
        change_point = 0

        df["before_onset"] = np.maximum(0, change_point - df[years_since_onset_col])
        df["after_onset"] = np.maximum(0, df[years_since_onset_col] - change_point)

        formula = f"{protein} ~ {formula_rhs}"
        model = smf.mixedlm(formula, df, groups=df[subject_id_col])
        result = model.fit()

        df_pred = df.copy()
        for cov in covariates:
            df_pred[cov] = mean_covariates[cov]
        df_pred["global_prediction"] = result.predict(df_pred)

        global_trend = df_pred.groupby(years_since_onset_col)["global_prediction"].mean()

        plt.plot(global_trend.index,
                 global_trend.values,
                 label=group_labels[group_key],
                 color=group_colors[group_key],
                 linewidth=2.5)

    # Plotting
    plt.figure(figsize=(6, 4), dpi=300)

    # Fit models
    fit_and_plot(df_status_change,
                 "before_onset + after_onset + " + " + ".join(covariates),
                 "status_change")

    fit_and_plot(df_normal,
                 "before_onset + " + " + ".join(covariates),
                 "normal")

    fit_and_plot(df_abnormal,
                 "after_onset + " + " + ".join(covariates),
                 "abnormal")

    plt.axvline(x=0, color="black", linestyle="--", linewidth=1.5, label="Onset")

    plt.xlabel("Years Since Onset")
    plt.ylabel(f"{protein} Expression")
    plt.title(f"{protein} Trajectory Aligned to Onset")

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(1)
        spine.set_color("black")

    # Dynamic ticks
    all_x = pd.concat([
        df_status_change[years_since_onset_col],
        df_normal[years_since_onset_col],
        df_abnormal[years_since_onset_col]
    ])
    min_x = int(np.floor(all_x.min() / 5) * 5)
    max_x = int(np.ceil(all_x.max() / 5) * 5)
    plt.xticks(np.arange(min_x, max_x + 1, 5))

    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True)
    plt.tight_layout()

    if export:
        os.makedirs(export_dir, exist_ok=True)
        for fmt in export_formats:
            plt.savefig(os.path.join(export_dir, f"{protein}_cplmm.{fmt}"),
                        format=fmt, bbox_inches="tight")

    plt.show()
