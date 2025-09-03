import pandas as pd

def compute_event_df(df, threshold, biomarker, group_label, onset_source='ONSET_AGE'):
    """
    Compute event data for Kaplan-Meier survival analysis based on biomarker threshold crossing.

    Parameters
    ----------
    df : pd.DataFrame
        Patient-level dataframe with longitudinal visits.
    threshold : float
        Biomarker threshold for defining event occurrence.
    biomarker : str
        Biomarker column name.
    group_label : str
        Label for patient group (e.g., 'Normal', 'MCI').
    onset_source : str
        Column name for onset age reference (default: 'ONSET_AGE').

    Returns
    -------
    pd.DataFrame
        Event dataframe with columns: SUBID, time (years since onset), event (0/1), group.
    """
    event_data = []
    for subid, g in df.groupby("SUBID"):
        g = g.sort_values("PROCEDURE_AGE")
        onset_age = g[onset_source].iloc[0]
        crossed = g[g[biomarker] >= threshold]
        if not crossed.empty:
            time = crossed["PROCEDURE_AGE"].iloc[0] - onset_age
            event = 1
        else:
            time = g["PROCEDURE_AGE"].max() - onset_age
            event = 0
        event_data.append({"SUBID": subid, "time": time, "event": event, "group": group_label})
    return pd.DataFrame(event_data)
