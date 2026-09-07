"""Load cell-count.csv into a normalized SQLite database.

Usage:
    python load_data.py [csv_path] [db_path]

With no arguments, reads ./cell-count.csv and writes ./cell_counts.db
(both resolved relative to this script's directory, i.e. the repo root).

Re-running is safe: the schema is dropped and rebuilt inside a single
transaction, so a failed run leaves the previous database untouched and a
successful run always reflects exactly the current CSV contents.
"""

import sqlite3
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_CSV_PATH = REPO_ROOT / "cell-count.csv"
DEFAULT_DB_PATH = REPO_ROOT / "cell_counts.db"

# Cell population columns in the source CSV, in the order they should be
# loaded into the `populations` lookup table.
CELL_POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

REQUIRED_COLUMNS = [
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    *CELL_POPULATIONS,
]

# Columns that must be present and non-blank on every row.
NON_EMPTY_COLUMNS = ["project", "subject", "condition", "sex", "treatment", "sample", "sample_type"]

# Columns that must parse as non-negative integers.
INTEGER_COLUMNS = ["age", "time_from_treatment_start", *CELL_POPULATIONS]

SUBJECT_ATTRIBUTE_COLUMNS = ["project", "condition", "sex", "age", "treatment", "response"]

SCHEMA_SQL = """
DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS populations;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id   TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id   TEXT PRIMARY KEY,
    project_id   TEXT NOT NULL REFERENCES projects(project_id),
    condition    TEXT NOT NULL,
    sex          TEXT NOT NULL,
    age          INTEGER NOT NULL,
    treatment    TEXT,
    response     TEXT
);

CREATE TABLE samples (
    sample_id                  TEXT PRIMARY KEY,
    subject_id                 TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type                TEXT NOT NULL,
    time_from_treatment_start  INTEGER NOT NULL
);

CREATE TABLE populations (
    population_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    name            TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id       TEXT NOT NULL REFERENCES samples(sample_id),
    population_id   INTEGER NOT NULL REFERENCES populations(population_id),
    count           INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population_id)
);

CREATE INDEX idx_subjects_project ON subjects(project_id);
CREATE INDEX idx_samples_subject ON samples(subject_id);
CREATE INDEX idx_cell_counts_population ON cell_counts(population_id);
"""


class MalformedDataError(ValueError):
    """Raised when the CSV is missing, unparsable, or fails validation."""


def _line_numbers(df, mask):
    """Map a boolean mask over `df` to 1-indexed CSV line numbers (header is line 1)."""
    return df.index[mask].tolist()


def read_and_validate_csv(csv_path):
    """Read and validate the CSV with pandas. Returns a cleaned DataFrame
    indexed by source CSV line number, with `treatment`/`response` holding
    Python/NaN nulls where the source said "none" or left the field blank,
    and all INTEGER_COLUMNS cast to int. 
    """
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Could not find CSV file at {csv_path}. "
            f"Expected cell-count.csv in the repository root."
        )

    try:
        df = pd.read_csv(csv_path, dtype=str)
    except pd.errors.EmptyDataError:
        raise MalformedDataError(f"{csv_path} is empty (no header row)") from None
    except pd.errors.ParserError as exc:
        raise MalformedDataError(f"{csv_path} could not be parsed as CSV: {exc}") from None

    # Header is CSV line 1, so the first data row is line 2.
    df.index = range(2, len(df) + 2)

    missing_columns = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_columns:
        raise MalformedDataError(
            f"{csv_path} is missing column(s) {missing_columns} -- check the header row"
        )
    df = df[REQUIRED_COLUMNS].copy()

    for col in NON_EMPTY_COLUMNS:
        blank_mask = df[col].isna() | (df[col].str.strip() == "")
        if blank_mask.any():
            raise MalformedDataError(
                f"cell-count.csv line(s) {_line_numbers(df, blank_mask)}: "
                f"missing required field '{col}'"
            )
        df[col] = df[col].str.strip()

    for col in INTEGER_COLUMNS:
        numeric = pd.to_numeric(df[col], errors="coerce")
        bad_mask = numeric.isna() | (numeric != numeric.round())
        if bad_mask.any():
            raise MalformedDataError(
                f"cell-count.csv line(s) {_line_numbers(df, bad_mask)}: "
                f"field '{col}' must be an integer"
            )
        negative_mask = numeric < 0
        if negative_mask.any():
            raise MalformedDataError(
                f"cell-count.csv line(s) {_line_numbers(df, negative_mask)}: "
                f"field '{col}' must be >= 0"
            )
        df[col] = numeric.astype(int)

    # "none" treatment and blank response represent absence, not text.
    df["treatment"] = df["treatment"].where(df["treatment"].str.lower() != "none", None)
    df["response"] = df["response"].where(df["response"].notna(), None)

    return df


def create_schema(conn):
    conn.executescript(SCHEMA_SQL)


def load_dataframe(conn, df):
    """Populate all tables from the validated DataFrame via DataFrame.to_sql.
    Raises MalformedDataError on inconsistent data (e.g. the same subject
    appearing with conflicting attributes, or a duplicate sample id).
    """
    projects_df = pd.DataFrame({"project_id": df["project"].unique()})
    projects_df.to_sql("projects", conn, if_exists="append", index=False)

    conflicting = df.groupby("subject")[SUBJECT_ATTRIBUTE_COLUMNS].nunique(dropna=False)
    conflicting = conflicting.index[(conflicting > 1).any(axis=1)].tolist()
    if conflicting:
        raise MalformedDataError(
            f"cell-count.csv: subject(s) {conflicting} have conflicting "
            f"metadata across rows (expected immutable per subject)"
        )

    subjects_df = (
        df[["subject", *SUBJECT_ATTRIBUTE_COLUMNS]]
        .drop_duplicates(subset="subject")
        .rename(columns={"subject": "subject_id", "project": "project_id"})
    )
    subjects_df.to_sql("subjects", conn, if_exists="append", index=False)

    dup_mask = df["sample"].duplicated()
    if dup_mask.any():
        raise MalformedDataError(
            f"cell-count.csv line(s) {_line_numbers(df, dup_mask)}: "
            f"duplicate sample id(s) {sorted(set(df['sample'][dup_mask]))}"
        )

    samples_df = (
        df[["sample", "subject", "sample_type", "time_from_treatment_start"]]
        .rename(columns={"sample": "sample_id", "subject": "subject_id"})
    )
    samples_df.to_sql("samples", conn, if_exists="append", index=False)

    populations_df = pd.DataFrame({"name": CELL_POPULATIONS})
    populations_df.to_sql("populations", conn, if_exists="append", index=False)
    population_ids = pd.read_sql("SELECT population_id, name FROM populations", conn)

    cell_counts_df = (
        df[["sample", *CELL_POPULATIONS]]
        .melt(id_vars="sample", var_name="name", value_name="count")
        .rename(columns={"sample": "sample_id"})
        .merge(population_ids, on="name")[["sample_id", "population_id", "count"]]
    )
    cell_counts_df.to_sql("cell_counts", conn, if_exists="append", index=False)


def print_summary(conn, db_path):
    tables = ["projects", "subjects", "samples", "populations", "cell_counts"]
    counts = {t: conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0] for t in tables}
    print(f"Loaded cell-count.csv into {db_path}")
    for table, count in counts.items():
        print(f"  {table}: {count}")


def main():
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_CSV_PATH
    db_path = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_DB_PATH

    try:
        df = read_and_validate_csv(csv_path)
    except (FileNotFoundError, MalformedDataError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        conn.execute("BEGIN")
        create_schema(conn)
        load_dataframe(conn, df)
        conn.commit()
    except Exception as exc:
        conn.rollback()
        print(f"Error loading data, database left unchanged: {exc}", file=sys.stderr)
        sys.exit(1)
    else:
        print_summary(conn, db_path)
    finally:
        conn.close()


if __name__ == "__main__":
    main()
