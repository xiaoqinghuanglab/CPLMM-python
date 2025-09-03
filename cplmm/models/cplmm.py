import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from sklearn.metrics import mean_squared_error
from typing import List, Optional

def fit_cplmm_all_proteins(
    df_status_change: pd.DataFrame,
    df_normal: pd.DataFrame,
    df_abnormal: pd.DataFrame,
    protein_list: List[str],
    covariates: List[str] = ["SEX", "BASELINE_AGE"],
    subject_id_col: str = "SUBID",
    years_since_onset_col: str = "years_since_onset"
) -> pd.DataFrame:
    """
    Fit change-point LMM across multiple proteins for status-change, normal, and abnormal groups.
    
    Returns a dataframe with slopes, SEs, and model metrics.
    """

    results = {}

    def fit_lmm_and_get_slopes(df, formula, label):
        if df.empty:
            return np.nan, np.nan, None

        df = df.copy()
        cp = 0
        df['before_onset'] = np.maximum(0, cp - df[years_since_onset_col])
        df['after_onset'] = np.maximum(0, df[years_since_onset_col] - cp)

        model = smf.mixedlm(formula, df, groups=df[subject_id_col])
        result = model.fit()

        fe = result.fe_params
        bse = result.bse

        slope_before = -fe.get("before_onset", np.nan)
        slope_after = fe.get("after_onset", np.nan)
        se_before = bse.get("before_onset", np.nan)
        se_after = bse.get("after_onset", np.nan)

        return slope_before, slope_after, result, se_before, se_after

    for protein in protein_list:
        formula_status = f"{protein} ~ before_onset + after_onset + {' + '.join(covariates)}"
        formula_normal = f"{protein} ~ before_onset + {' + '.join(covariates)}"
        formula_abnormal = f"{protein} ~ after_onset + {' + '.join(covariates)}"

        s1, s3, r_status, se1, se3 = fit_lmm_and_get_slopes(df_status_change, formula_status, "Status Change")
        s2, _, r_normal, se2, _ = fit_lmm_and_get_slopes(df_normal, formula_normal, "Normal")
        _, s4, r_abnormal, _, se4 = fit_lmm_and_get_slopes(df_abnormal, formula_abnormal, "Abnormal")

        def get_metrics(result, df, group_name):
            if result is None:
                return np.nan, np.nan, np.nan
            aic = result.aic
            bic = result.bic
            y_true = df[protein]
            y_pred = result.predict(df)
            mse = mean_squared_error(y_true, y_pred)
            return aic, bic, mse

        aic_s, bic_s, mse_s = get_metrics(r_status, df_status_change, "status")
        aic_n, bic_n, mse_n = get_metrics(r_normal, df_normal, "normal")
        aic_a, bic_a, mse_a = get_metrics(r_abnormal, df_abnormal, "abnormal")

        results[protein] = {
            'Protein': protein,
            'Beta 1': s1, 'SE Beta 1': se1,
            'Beta 2': s2, 'SE Beta 2': se2,
            'Beta 3': s3, 'SE Beta 3': se3,
            'Beta 4': s4, 'SE Beta 4': se4,
            'Intercept': r_status.fe_params.get("Intercept", np.nan) if r_status else np.nan,
            'AIC Status': aic_s, 'BIC Status': bic_s, 'MSE Status': mse_s,
            'AIC Normal': aic_n, 'BIC Normal': bic_n, 'MSE Normal': mse_n,
            'AIC Abnormal': aic_a, 'BIC Abnormal': bic_a, 'MSE Abnormal': mse_a
        }

    df_results = pd.DataFrame.from_dict(results, orient='index')
    df_results.reset_index(drop=True, inplace=True)
    return df_results
