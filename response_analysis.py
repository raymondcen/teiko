"""Compare immune cell population relative frequencies between miraclib
responders and non-responders, restricted to melanoma / PBMC samples.

Usage:
    python response_analysis.py [db_path]

Cohort: melanoma subjects on miraclib, PBMC samples, response recorded,
restricted to the time_from_treatment_start = 0 (baseline/pre-treatment)
draw. Two reasons for the baseline restriction rather than pooling all
three draws (t=0/7/14) per subject:

  1. Predictive vs. reactive: to *predict* response from a cell population,
     the measurement has to precede the outcome. Using post-dose samples
     (t=7, t=14) risks measuring an effect of the treatment/response
     itself rather than a baseline biomarker of it.
  2. Independence: each subject has exactly one PBMC sample at t=0, so
     each subject contributes exactly one observation. Pooling all three
     timepoints would pseudoreplicate (3 correlated samples per subject
     treated as independent), inflating apparent significance.

Prints a per-population statistics table (Mann-Whitney U test, two-sided,
with Benjamini-Hochberg FDR correction across the five populations since
five hypotheses are tested at once) and writes an interactive boxplot to
response_boxplots.html.
"""

import json
import sys
from pathlib import Path

import pandas as pd
import plotly.express as px
import sqlite3
from scipy.stats import false_discovery_control, mannwhitneyu

from analysis import relative_frequencies

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_DB_PATH = REPO_ROOT / "cell_counts.db"
DEFAULT_PLOT_PATH = REPO_ROOT / "response_boxplots.html"

# Cohort definition: melanoma subjects on miraclib, PBMC samples only, with
# a recorded response, restricted to the baseline (pre-treatment) draw --
# see the module docstring for why.
COHORT_SQL = """
SELECT
    sm.sample_id AS sample,
    sj.response AS response
FROM samples sm
JOIN subjects sj ON sj.subject_id = sm.subject_id
WHERE sj.condition = 'melanoma'
  AND sj.treatment = 'miraclib'
  AND sm.sample_type = 'PBMC'
  AND sm.time_from_treatment_start = 0
  AND sj.response IN ('yes', 'no');
"""


def cohort_frequencies(conn):
    """Relative frequency rows (sample, total_count, population, count,
    percentage) restricted to melanoma/miraclib/PBMC samples with a
    recorded response, plus a `response` column ("yes"/"no").
    """
    cohort = pd.read_sql(COHORT_SQL, conn)
    freq = relative_frequencies(conn)
    return freq.merge(cohort, on="sample", how="inner")


def compare_responders(df):
    """One row per population: group sizes, medians, Mann-Whitney U
    statistic, raw p-value, and BH-adjusted p-value across all populations.
    """
    rows = []
    for population, group in df.groupby("population"):
        yes = group.loc[group["response"] == "yes", "percentage"]
        no = group.loc[group["response"] == "no", "percentage"]
        stat, p_value = mannwhitneyu(yes, no, alternative="two-sided")
        # Rank-biserial correlation: effect size on a [-1, 1] scale, robust
        # to the same non-normality Mann-Whitney tolerates.
        effect_size = 1 - (2 * stat) / (len(yes) * len(no))
        rows.append(
            {
                "population": population,
                "n_responders": len(yes),
                "n_non_responders": len(no),
                "median_responders_pct": yes.median(),
                "median_non_responders_pct": no.median(),
                "effect_size": effect_size,
                "u_statistic": stat,
                "p_value": p_value,
            }
        )
    result = pd.DataFrame(rows).sort_values("population").reset_index(drop=True)
    result["p_adj"] = false_discovery_control(result["p_value"], method="bh")
    result["significant"] = result["p_adj"] < 0.05
    return result


def build_boxplot(df):
    ordered = sorted(df["population"].unique())
    fig = px.box(
        df,
        x="response",
        y="percentage",
        color="response",
        facet_col="population",
        category_orders={"response": ["no", "yes"], "population": ordered},
        points="outliers",
        labels={
            "percentage": "Relative frequency (%)",
            "response": "Response to miraclib",
        },
        title="Melanoma / PBMC: cell population relative frequency by miraclib response",
    )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.update_layout(showlegend=False)
    return fig


def write_boxplot(df, plot_path):
    build_boxplot(df).write_html(plot_path, include_plotlyjs="cdn")


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
        df = cohort_frequencies(conn)
    finally:
        conn.close()

    stats = compare_responders(df)
    print(stats.to_string(index=False))

    write_boxplot(df, DEFAULT_PLOT_PATH)
    print(f"\nWrote boxplot to {DEFAULT_PLOT_PATH}")


if __name__ == "__main__":
    main()
