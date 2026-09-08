"""Summarize the baseline miraclib/melanoma/PBMC cohort.

Usage:
    python baseline_cohort.py [db_path]

Identifies all melanoma, PBMC samples at baseline
(time_from_treatment_start = 0) from subjects treated with miraclib, then
reports three breakdowns: samples per project, subjects per response, and
subjects per sex.
"""

import sqlite3
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_DB_PATH = REPO_ROOT / "cell_counts.db"

BASELINE_COHORT_SQL = """
SELECT
    sm.sample_id AS sample,
    sj.subject_id AS subject,
    sj.project_id AS project,
    sj.sex AS sex,
    sj.response AS response
FROM samples sm
JOIN subjects sj ON sj.subject_id = sm.subject_id
WHERE sj.condition = 'melanoma'
  AND sj.treatment = 'miraclib'
  AND sm.sample_type = 'PBMC'
  AND sm.time_from_treatment_start = 0;
"""


def baseline_cohort(conn):
    """One row per baseline sample in the melanoma/miraclib/PBMC cohort,
    with the owning subject's project, sex, and response.
    """
    return pd.read_sql(BASELINE_COHORT_SQL, conn)


def samples_per_project(df):
    return (
        df.groupby("project")["sample"]
        .nunique()
        .rename("n_samples")
        .reset_index()
        .sort_values("project")
    )


def subjects_per_response(df):
    subjects = df.drop_duplicates(subset="subject")
    return (
        subjects["response"]
        .value_counts(dropna=False)
        .rename_axis("response")
        .reset_index(name="n_subjects")
    )


def subjects_per_sex(df):
    subjects = df.drop_duplicates(subset="subject")
    return (
        subjects["sex"]
        .value_counts(dropna=False)
        .rename_axis("sex")
        .reset_index(name="n_subjects")
    )


def main():
    db_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_DB_PATH

    if not db_path.exists():
        print(
            f"Error: could not find database at {db_path}. "
            f"Run `python load_data.py` first.",
            file=sys.stderr,
        )
        sys.exit(1)

    conn = sqlite3.connect(db_path)
    try:
        df = baseline_cohort(conn)
    finally:
        conn.close()

    n_samples = df["sample"].nunique()
    n_subjects = df["subject"].nunique()
    print(f"Baseline melanoma/miraclib/PBMC cohort: {n_samples} samples, {n_subjects} subjects\n")

    print("Samples per project:")
    print(samples_per_project(df).to_string(index=False))

    print("\nSubjects per response:")
    print(subjects_per_response(df).to_string(index=False))

    print("\nSubjects per sex:")
    print(subjects_per_sex(df).to_string(index=False))


if __name__ == "__main__":
    main()
