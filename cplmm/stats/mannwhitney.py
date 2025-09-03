from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import pandas as pd

def compare_groups_mannwhitney(
    combined_expr: pd.DataFrame,
    gene_list: list,
    group1: str,
    group2: str,
    alpha: float = 0.05,
    rank_by: str = "FDR"
) -> pd.DataFrame:
    """
    Perform Mann-Whitney U test to compare expression between two groups for a set of genes.

    Parameters
    ----------
    combined_expr : pd.DataFrame
        Long-format expression dataframe (columns: Gene, Expression, Source).
    gene_list : list
        List of genes/proteins to test.
    group1 : str
        First group label in `Source` column.
    group2 : str
        Second group label in `Source` column.
    alpha : float
        Significance threshold (for marking significance after FDR correction).
    rank_by : str
        Column to sort results by ("FDR", "p_value", or "Delta_mean").

    Returns
    -------
    pd.DataFrame
        Results dataframe with columns:
        - Gene
        - group means
        - Delta_mean
        - U_statistic
        - p_value
        - FDR (adjusted p-value)
        - Significant (boolean flag based on FDR < alpha)
    """
    results = []

    for gene in gene_list:
        expr1 = combined_expr[(combined_expr["Gene"] == gene) & (combined_expr["Source"] == group1)]["Expression"].dropna()
        expr2 = combined_expr[(combined_expr["Gene"] == gene) & (combined_expr["Source"] == group2)]["Expression"].dropna()

        if len(expr1) > 0 and len(expr2) > 0:
            stat, p = mannwhitneyu(expr1, expr2, alternative="two-sided")
            delta = expr2.mean() - expr1.mean()

            results.append({
                "Gene": gene,
                f"{group1}_mean": expr1.mean(),
                f"{group2}_mean": expr2.mean(),
                "Delta_mean": delta,
                "U_statistic": stat,
                "p_value": p
            })

    results_df = pd.DataFrame(results)

    # FDR correction
    results_df["FDR"] = multipletests(results_df["p_value"], method="fdr_bh")[1]
    results_df["Significant"] = results_df["FDR"] < alpha

    # Sort results
    if rank_by in results_df.columns:
        results_df = results_df.sort_values(rank_by, ascending=(rank_by != "Delta_mean")).reset_index(drop=True)
    else:
        results_df = results_df.sort_values("FDR").reset_index(drop=True)

    return results_df
