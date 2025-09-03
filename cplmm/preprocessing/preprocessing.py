import pandas as pd

def calculate_years_since_onset(df: pd.DataFrame,
                                 age_col: str = "age",
                                 onset_age_col: str = "onset_age",
                                 new_col: str = "years_since_onset") -> pd.DataFrame:
    """
    Calculate time since disease onset.
    """
    df[new_col] = df[age_col] - df[onset_age_col]
    return df


def add_piecewise_age(df: pd.DataFrame,
                      age_col: str = "age",
                      onset_age_col: str = "onset_age",
                      new_col: str = "piecewise_age") -> pd.DataFrame:
    """
    Create piecewise age variable using onset age as the knot.
    """
    df[new_col] = df[age_col] - df[onset_age_col]
    return df


def set_onset_age_for_normals(df: pd.DataFrame,
                              subject_id_col: str = "subject_id",
                              status_col: str = "status_raw",
                              age_col: str = "age",
                              mutated_col: str = "decage_mutated") -> pd.DataFrame:
    """
    For patients who remain 'Normal' across visits, set onset age to max age.
    """
    normal_ids = (
        df.groupby(subject_id_col)
          .filter(lambda x: set(x[status_col]) == {'Normal'})[subject_id_col]
          .unique()
    )

    for subid in normal_ids:
        max_age = df.loc[df[subject_id_col] == subid, age_col].max()
        df.loc[df[subject_id_col] == subid, mutated_col] = max_age

    return df


def correct_status_by_onset(df: pd.DataFrame,
                            age_col: str = "age",
                            onset_age_col: str = "onset_age",
                            status_col: str = "status_cleaned") -> pd.DataFrame:
    """
    Fix logical inconsistencies based on onset age.
    """
    df.loc[(df[onset_age_col] > df[age_col]) & (df[status_col] == 'Abnormal'), status_col] = 'Normal'
    df.loc[(df[onset_age_col] < df[age_col]) & (df[status_col] == 'Normal'), status_col] = 'Abnormal'
    return df


def enforce_unidirectional_status_change(df: pd.DataFrame,
                                         subject_id_col: str = "subject_id",
                                         status_col: str = "status_cleaned",
                                         date_col: str = "procedure_date") -> pd.DataFrame:
    """
    Once abnormal, always abnormal.
    """
    df = df.sort_values(by=[subject_id_col, date_col])

    def correct_status(sub_df):
        abnormal_seen = False
        for idx, row in sub_df.iterrows():
            if row[status_col] == "Abnormal":
                abnormal_seen = True
            if abnormal_seen:
                sub_df.at[idx, status_col] = "Abnormal"
        return sub_df

    return df.groupby(subject_id_col, group_keys=False).apply(correct_status)


def identify_status_change_subjects(df: pd.DataFrame,
                                    subject_id_col: str = "subject_id",
                                    status_col: str = "status_cleaned",
                                    date_col: str = "procedure_date") -> pd.DataFrame:
    """
    Return only subjects that changed from Normal to Abnormal.
    """
    df_sorted = df.sort_values(by=[subject_id_col, date_col])
    return df_sorted.groupby(subject_id_col).filter(
        lambda x: (x[status_col] == "Normal").any() and (x[status_col] != "Normal").any()
    )


def get_status_groups(df: pd.DataFrame,
                      subject_id_col: str = "subject_id",
                      status_col: str = "status_cleaned") -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Return normal-only and abnormal-only groups.
    """
    normal_ids = df.groupby(subject_id_col)[status_col].apply(lambda x: (x == "Normal").all())
    abnormal_ids = df.groupby(subject_id_col)[status_col].apply(lambda x: (x != "Normal").all())

    df_normal = df[df[subject_id_col].isin(normal_ids[normal_ids].index)]
    df_abnormal = df[df[subject_id_col].isin(abnormal_ids[abnormal_ids].index)]

    return df_normal, df_abnormal
