"""Interactive dashboard presenting the Part 2-4 analysis results.

Usage:
    streamlit run dashboard.py

Reads cell_counts.db (built by `python load_data.py`, i.e. `make
pipeline`) and computes every table/chart live from it via the same
functions used by analysis.py, response_analysis.py, and
baseline_cohort.py -- no precomputed/cached results, so the dashboard
always reflects the current database.
"""

import sqlite3
from pathlib import Path

import streamlit as st

from analysis import relative_frequencies
from baseline_cohort import (
    baseline_cohort,
    samples_per_project,
    subjects_per_response,
    subjects_per_sex,
)
from response_analysis import build_boxplot, compare_responders, cohort_frequencies

REPO_ROOT = Path(__file__).resolve().parent
DB_PATH = REPO_ROOT / "cell_counts.db"

st.set_page_config(page_title="Teiko cell-count analysis", layout="wide")
st.title("Teiko cell-count analysis")

if not DB_PATH.exists():
    st.error(f"No database found at `{DB_PATH}`. Run `make pipeline` first.")
    st.stop()

conn = sqlite3.connect(DB_PATH)

st.header("Part 2 -- Relative frequency by sample")
st.caption(
    "One row per (sample, population): total cell count, count, and "
    "percentage of the sample's total."
)
freq_df = relative_frequencies(conn)
samples = st.multiselect(
    "Filter to sample(s) (leave empty to show all)",
    options=sorted(freq_df["sample"].unique()),
)
st.dataframe(
    freq_df[freq_df["sample"].isin(samples)] if samples else freq_df,
    width="stretch",
    hide_index=True,
)

st.header("Part 3 -- Miraclib response signal (melanoma, PBMC, baseline)")
st.caption(
    "Baseline (time_from_treatment_start = 0) PBMC samples from melanoma "
    "subjects on miraclib, compared responders vs. non-responders per "
    "population with a two-sided Mann-Whitney U test and Benjamini-Hochberg "
    "FDR correction across the five populations."
)
response_df = cohort_frequencies(conn)
stats_df = compare_responders(response_df)
st.plotly_chart(build_boxplot(response_df), width="stretch")
st.dataframe(stats_df, width="stretch", hide_index=True)
if not stats_df["significant"].any():
    st.info(
        "No population reaches significance after FDR correction "
        f"(lowest q = {stats_df['p_adj'].min():.3f})."
    )

st.header("Part 4 -- Baseline miraclib / melanoma / PBMC cohort")
cohort_df = baseline_cohort(conn)
col1, col2, col3 = st.columns(3)
with col1:
    st.subheader("Samples per project")
    st.dataframe(samples_per_project(cohort_df), width="stretch", hide_index=True)
with col2:
    st.subheader("Subjects per response")
    st.dataframe(subjects_per_response(cohort_df), width="stretch", hide_index=True)
with col3:
    st.subheader("Subjects per sex")
    st.dataframe(subjects_per_sex(cohort_df), width="stretch", hide_index=True)

conn.close()
