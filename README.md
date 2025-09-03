# 🧬 Change Point Linear MIxed Effect Model (CPLMM)

**cplmm** is a Python package for longitudinal proteomics and biomarker analysis, designed for research involving disease progression, biomarker discovery, and pathway enrichment.
It provides tools for preprocessing, change-point modeling (CPLMM), statistical testing, survival analysis, and publication-ready visualizations.

---

## ✨ Features

- 📊 **Preprocessing**: Automated onset age calculation, patient filtering, and data preparation
- 📈 **Change-point Linear Mixed Models (CPLMM)**: Detect biomarker trajectory shifts
- 🧪 **Statistical Tests**: Wald tests (with FDR correction) and Mann-Whitney U comparisons
- ⏳ **Survival Analysis**: Kaplan-Meier plots, at-risk tables, and log-rank tests
- 🔬 **Pathway Analysis Visualization**: Bubble plots, gene-pathway heatmaps, and top pathways
- 🎨 **Publication-ready Figures**: Plots with customizable palettes and fonts

---

## 🔧 Installation

Clone the repository and install locally:

```bash
git clone https://github.com/xiaoqinghuanglab/CPLMM-python.git
cd cplmm
pip install -e .
```

---

## 📂 Package Structure

```
cplmm/
├── preprocessing/          # Data cleaning & expression prep
├── models/                 # Change-point LMM (CPLMM)
├── stats/                  # Wald tests & Mann-Whitney comparisons
├── survival/               # Kaplan-Meier and event computation
├── viz/                    # Visualization (CPLMM plots, volcano, pathways, etc.)
├── utils/                  # JAMA palettes & styling
└── __init__.py             # Unified import access
```

## 🗂 Required DataFrame Structure

The package expects a **longitudinal biomarker dataframe** with standardized patient-level and biomarker-level data.  
For flexibility, all functions allow mapping your dataset’s custom column names to general arguments.

### 🔑 **Core Columns (Required)**
| **General Argument**     | **Description**                                                                 |
|--------------------------|---------------------------------------------------------------------------------|
| `subject_id`            | Unique patient identifier (string or integer).                                  |
| `age_col`               | Patient's age at sample collection (float).                                     |
| `date_col`              | Date of sample collection (datetime).                                           |
| `status_col_raw`        | Initial diagnostic status (e.g., `Normal`, `Abnormal`).                         |
| `status_col_cleaned`    | Unidirectional status corrected for onset progression (derived field). [Can be created using the Pre Processing functions of the package.]         |
| `onset_age_col`         | Age of disease onset (calculated or provided).                                  |
| `decage_mutated_col`    | Adjusted diagnosis age for patients remaining Normal throughout follow-up.      |

All **protein/biomarker expression values** should be in separate columns, with numeric values representing measured expression.

---

## 🛠 Preprocessing Utilities

`cplmm.preprocessing` provides tools to prepare longitudinal biomarker datasets for analysis.  
These functions handle **time alignment, onset-based transformations, status corrections, and group filtering**.

### 🔑 **Available Functions**

1. **`calculate_years_since_onset`**  
   - Computes years relative to disease onset for each visit.  
   - Adds a new column (default: `years_since_onset`).

2. **`add_piecewise_age`**  
   - Creates a piecewise variable (age relative to onset) for changepoint modeling.

3. **`set_onset_age_for_normals`**  
   - For patients who remain "Normal" throughout follow-up, assigns onset age as their maximum recorded age.

4. **`correct_status_by_onset`**  
   - Ensures logical consistency between diagnosis status and onset age.

5. **`enforce_unidirectional_status_change`**  
   - Ensures that once a subject transitions to "Abnormal," they remain abnormal in all subsequent visits.

6. **`identify_status_change_subjects`**  
   - Filters to patients who transition from "Normal" to "Abnormal."

7. **`get_status_groups`**  
   - Splits the dataset into **normal-only** and **abnormal-only** subjects.

---

### 🖥 **Example Usage**
```python
import pandas as pd
from cplmm.preprocessing import (
    calculate_years_since_onset,
    add_piecewise_age,
    set_onset_age_for_normals,
    correct_status_by_onset,
    enforce_unidirectional_status_change,
    identify_status_change_subjects,
    get_status_groups
)

# Load your dataframe
df = pd.read_csv("data.csv")

# 1. Calculate years since onset and piecewise age
df = calculate_years_since_onset(df, age_col="PROCEDURE_AGE", onset_age_col="ONSET_AGE")
df = add_piecewise_age(df, age_col="PROCEDURE_AGE", onset_age_col="ONSET_AGE")

# 2. Adjust onset age for normal-only patients
df = set_onset_age_for_normals(df, subject_id_col="SUBID", status_col="Status", age_col="PROCEDURE_AGE", mutated_col="DECAGE_Mutated")

# 3. Correct statuses and enforce unidirectional change
df = correct_status_by_onset(df, age_col="PROCEDURE_AGE", onset_age_col="ONSET_AGE", status_col="Status")
df = enforce_unidirectional_status_change(df, subject_id_col="SUBID", status_col="Status", date_col="PROCEDURE_DATE")

# 4. Identify status-change patients
status_change_df = identify_status_change_subjects(df, subject_id_col="SUBID", status_col="Status", date_col="PROCEDURE_DATE")

# 5. Get normal-only and abnormal-only groups
df_normal, df_abnormal = get_status_groups(df, subject_id_col="SUBID", status_col="Status")
```
---
## 🚀 Quick Start

### Preprocessing Example
```python
from cplmm import prepare_combined_expression

combined_expr = prepare_combined_expression(
    working_df, df_normal_only, df_abnormal_only,
    top_genes=top_100[:50]
)

### Change-Point Linear Mixed Model (CPLMM)

from cplmm import fit_cplmm

results_df = fit_cplmm(
    working_df, df_normal_only, df_abnormal_only,
    proteins=protein_list,
    covariates=["SEX", "PROCEDURE_AGE"]
)
```

### Wald Test

```python
from cplmm import wald_test
wald_df = wald_test(results_df)
```

### Mann-Whitney U Comparison

```python
from cplmm import compare_groups_mannwhitney

results_mci = compare_groups_mannwhitney(
    combined_expr, gene_list=top_100[:50],
    group1="Normal_only", group2="MCI"
)
```

### Kaplan-Meier Survival Analysis

```python
from cplmm import plot_km_with_threshold

plot_km_with_threshold(
    biomarker_name="NFL",
    threshold=threshold_dict["NFL"],
    wd_df=working_df,
    normal_df=df_normal_only,
    abnm_df=df_abnormal_only
)
```

---

## 📊 Visualization Examples

**CPLMM Trajectories:**
```python
from cplmm import plot_cplmm
plot_cplmm(working_df, df_normal_only, df_abnormal_only, protein="NFL")
```

**Volcano Plot (Wald Test):**
```python
from cplmm import plot_wald_volcano
plot_wald_volcano(wald_df, p_col="P-value_1", fdr_col="Adjusted P-value (FDR)_1")
```

**Pathway Bubble Plot:**
```python
from cplmm import plot_pathway_bubble
plot_pathway_bubble(cleaned_pathway_df)
```
