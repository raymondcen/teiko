"""Compute per-sample relative frequencies of each immune cell population.

Usage:
    python analysis.py [db_path]

With no arguments, reads ./cell_counts.db (resolved relative to this
script's directory, i.e. the repo root) and prints the relative frequency
table to stdout.
"""

import sqlite3
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_DB_PATH = REPO_ROOT / "cell_counts.db"

# One row per (sample, population): total_count is the sum of all five
# population counts for that sample, computed via a window function so it
# sits alongside each population's own count/percentage.
RELATIVE_FREQUENCY_SQL = """
SELECT
    s.sample_id AS sample,
    SUM(cc.count) OVER (PARTITION BY s.sample_id) AS total_count,
    p.name AS population,
    cc.count AS count,
    ROUND(
        100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY s.sample_id), 2
    ) AS percentage
FROM cell_counts cc
JOIN samples s ON s.sample_id = cc.sample_id
JOIN populations p ON p.population_id = cc.population_id
ORDER BY s.sample_id, p.population_id;
"""


def relative_frequencies(conn):
    """Return a DataFrame with columns sample, total_count, population,
    count, percentage -- one row per (sample, population) pair.
    """
    return pd.read_sql(RELATIVE_FREQUENCY_SQL, conn)


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
        df = relative_frequencies(conn)
    finally:
        conn.close()

    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
