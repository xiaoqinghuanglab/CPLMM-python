import pandas as pd

def prepare_combined_expression(
    df_normal_only: pd.DataFrame,
    df_status_change: pd.DataFrame,
    df_abnormal_only: pd.DataFrame,
    df_all: pd.DataFrame,
    subset_genes: list,
    category_col: str = "CATEGORY",
    categories: list = ["SCD", "MCI", "AD Dementia", "FTD Dementia"],
    normal_label: str = "Normal_only",
    status_label: str = "Status_change",
    abnormal_label: str = "Abnormal_only"
) -> pd.DataFrame:
    """
    Prepare a long-format expression dataframe for selected genes across patient categories.

    Parameters
    ----------
    df_normal_only : pd.DataFrame
        Expression dataframe for normal-only patients.
    df_status_change : pd.DataFrame
        Expression dataframe for status-change patients.
    df_abnormal_only : pd.DataFrame
        Expression dataframe for abnormal-only patients.
    df_all : pd.DataFrame
        Full dataframe including category labels.
    subset_genes : list
        List of genes/proteins to include.
    category_col : str
        Column indicating diagnosis category (e.g., 'CATEGORY').
    categories : list
        Diagnostic categories to subset (e.g., ['SCD', 'MCI', 'AD Dementia', 'FTD Dementia']).
    normal_label : str
        Label for normal-only patients in output.
    status_label : str
        Label for status-change patients in output.
    abnormal_label : str
        Label for abnormal-only patients in output.

    Returns
    -------
    pd.DataFrame
        Long-format dataframe with columns: ['Gene', 'Expression', 'Source'].
    """

    def melt_subset(df, label):
        melted = df[subset_genes].melt(var_name="Gene", value_name="Expression")
        melted["Expression"] = pd.to_numeric(melted["Expression"], errors="coerce")
        melted["Source"] = label
        return melted

    # Core groups
    melted_normal = melt_subset(df_normal_only, normal_label)
    melted_status = melt_subset(df_status_change, status_label)
    melted_abnormal = melt_subset(df_abnormal_only, abnormal_label)

    # Diagnostic categories
    category_melted = []
    for cat in categories:
        subset_df = df_all[df_all[category_col] == cat]
        melted = subset_df[subset_genes].melt(var_name="Gene", value_name="Expression")
        melted["Expression"] = pd.to_numeric(melted["Expression"], errors="coerce")
        melted["Source"] = cat.replace(" ", "_")
        category_melted.append(melted)

    # Combine all
    combined_expr = pd.concat(
        [melted_normal, melted_status, melted_abnormal] + category_melted,
        ignore_index=True
    )

    # Order genes
    combined_expr["Gene"] = pd.Categorical(combined_expr["Gene"], categories=subset_genes, ordered=True)

    return combined_expr
