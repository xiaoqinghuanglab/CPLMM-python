import pandas as pd
import numpy as np
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

def compute_wald_test(
    results_df: pd.DataFrame,
    adjust_p: bool = True,
    rank_by: int = 1,
    alpha: float = 0.05
) -> pd.DataFrame:
    """
    Compute Wald tests comparing pre- vs post-onset slopes and return ranked results.

    Parameters
    ----------
    results_df : pd.DataFrame
        Output from `fit_cplmm_all_proteins`, containing Beta slopes and SEs.
    adjust_p : bool
        Whether to apply FDR correction (Benjamini-Hochberg).
    rank_by : int
        Whether to rank genes by `p-value 1` (Status-change slopes) or `p-value 2` (Normal vs Abnormal slopes).
    alpha : float
        Significance level for marking significant results.

    Returns
    -------
    pd.DataFrame
        Wald test results including p-values, FDR-adjusted p-values, significance flags, and ranking.
    """
    wald_results = []

    for _, row in results_df.iterrows():
        # Extract betas and SEs
        beta_1, se_beta_1 = row.get('Beta 1', np.nan), row.get('SE Beta 1', np.nan)
        beta_3, se_beta_3 = row.get('Beta 3', np.nan), row.get('SE Beta 3', np.nan)
        beta_2, se_beta_2 = row.get('Beta 2', np.nan), row.get('SE Beta 2', np.nan)
        beta_4, se_beta_4 = row.get('Beta 4', np.nan), row.get('SE Beta 4', np.nan)

        # Skip if missing
        if any(pd.isna([beta_1, se_beta_1, beta_3, se_beta_3, beta_2, se_beta_2, beta_4, se_beta_4])):
            continue

        # Wald statistics
        wald_stat_1 = ((beta_1 - beta_3) ** 2) / (se_beta_1 ** 2 + se_beta_3 ** 2)
        wald_stat_2 = ((beta_2 - beta_4) ** 2) / (se_beta_2 ** 2 + se_beta_4 ** 2)

        # p-values
        p_val_1 = 1 - chi2.cdf(wald_stat_1, df=1)
        p_val_2 = 1 - chi2.cdf(wald_stat_2, df=1)

        wald_results.append({
            "Protein": row["Protein"],
            "Beta 1": beta_1, "SE Beta 1": se_beta_1,
            "Beta 3": beta_3, "SE Beta 3": se_beta_3,
            "Wald Statistic 1": wald_stat_1, "P-value 1": p_val_1,
            "Beta 2": beta_2, "SE Beta 2": se_beta_2,
            "Beta 4": beta_4, "SE Beta 4": se_beta_4,
            "Wald Statistic 2": wald_stat_2, "P-value 2": p_val_2
        })

    wald_df = pd.DataFrame(wald_results)

    # FDR adjustment
    if adjust_p:
        for i in [1, 2]:
            rejected, p_adj, _, _ = multipletests(wald_df[f"P-value {i}"], method="fdr_bh")
            wald_df[f"Adjusted P-value {i}"] = p_adj
            wald_df[f"Significant {i}"] = p_adj < alpha

    # Rank by chosen p-value (adjusted if available)
    rank_col = f"Adjusted P-value {rank_by}" if adjust_p else f"P-value {rank_by}"
    wald_df = wald_df.sort_values(by=rank_col).reset_index(drop=True)
    wald_df["Rank"] = range(1, len(wald_df) + 1)

    return wald_df
