# teiko

Loads immune cell population counts (`cell-count.csv`) into a normalized
SQLite database, computes relative-frequency and miraclib-response
statistics from it, and serves the results as an interactive dashboard.

**Dashboard:** run locally with `make dashboard` — see [Dashboard](#dashboard).

## Quickstart (GitHub Codespaces)

Works out of the box in a GitHub Codespace (Python 3.13 is preinstalled) or
any local Python 3.13+ environment with `make`.

```bash
# (Codespaces: open the repo, then in the terminal)
make setup       # pip install -r requirements.txt
make pipeline    # rebuilds cell_counts.db from cell-count.csv
make dashboard   # starts the Streamlit dashboard on port 8501
```

`make dashboard` prints a `Local URL`/`Network URL`. In Codespaces, a
"port forwarded" notification appears (or check the **Ports** tab) — open
that URL to view the dashboard in the browser; Codespaces proxies it over
HTTPS automatically.

`requirements.txt` installs `pandas` (CSV loading), `scipy` (statistics),
`plotly` (charts), and `streamlit` (the dashboard) — one `pip install` for
the whole project.

## Reproducing each part from the command line

The dashboard computes everything live, but each part also runs standalone
and prints its own table, useful for checking a number without opening a
browser:

```bash
python load_data.py             # Part 1: cell-count.csv -> cell_counts.db
python analysis.py               # Part 2: relative frequency per sample/population
python response_analysis.py     # Part 3: responder vs. non-responder stats + response_boxplots.html
python baseline_cohort.py       # Part 4: baseline miraclib/melanoma/PBMC cohort breakdown
```

All four accept an optional `db_path` argument (and `load_data.py` a
`csv_path` too), resolved relative to the repo root by default. Re-running
`load_data.py` is safe — it wipes and rebuilds the schema inside a single
transaction, so the database always reflects exactly the current CSV, and
a failed run (bad CSV, missing file) leaves the previous database
untouched.

Inspect the database directly with the `sqlite3` CLI (preinstalled in
Codespaces):

```bash
sqlite3 cell_counts.db ".schema"
sqlite3 -header -column cell_counts.db "SELECT * FROM samples LIMIT 5;"
```

## Schema

Five tables, normalized to 3NF:

```
projects                subjects                     samples                        cell_counts              populations
---------                --------                     -------                        -----------              -----------
project_id (PK)  <----+  subject_id (PK)      <----+  sample_id (PK)          <----+  sample_id (PK, FK)  +--> population_id (PK)
                       |  project_id (FK)            |  subject_id (FK)             |  population_id (PK, FK)  name
                       +  condition                  +  sample_type                +  count
                          sex                            time_from_treatment_start
                          age
                          treatment
                          response
```

- **projects** — one row per project (`prj1`, `prj2`, ...).
- **subjects** — one row per subject/patient: demographics (`condition`,
  `sex`, `age`) and clinical assignment (`treatment`, `response`). These are
  treated as immutable per subject — the source data confirms every row for
  a given `subject` carries identical values for these fields.
  `treatment = NULL` and `response = NULL` represent subjects with no
  treatment (the CSV's `none`/blank values), rather than the string `"none"`,
  so absence is queryable with standard `IS NULL` semantics instead of a
  magic string.
- **samples** — one row per physical sample draw: which subject it came
  from, `sample_type` (PBMC/WB), and `time_from_treatment_start`. This is
  the mutable, repeatable part of the data — a subject has many samples
  over time.
- **populations** — a small lookup table (`b_cell`, `cd8_t_cell`,
  `cd4_t_cell`, `nk_cell`, `monocyte`) giving each cell population a stable
  `population_id`.
- **cell_counts** — the measurement fact table, in long/melted form: one
  row per `(sample_id, population_id)` pair holding a `count`. This is the
  key normalization decision — see rationale below.

Foreign keys: `subjects.project_id -> projects`,
`samples.subject_id -> subjects`,
`cell_counts.sample_id -> samples`,
`cell_counts.population_id -> populations`.
`PRAGMA foreign_keys = ON` is enabled at load time, so referential
integrity is enforced by SQLite itself.

## Design rationale

**Separating subjects from samples** mirrors the real-world entities: a
subject is measured once for demographics/treatment/outcome but sampled
repeatedly over time. Storing subject attributes once (rather than
repeated on every sample row, as in the raw CSV) removes update anomalies
and duplication, and makes "how many subjects" vs. "how many samples"
trivial, correct questions instead of `COUNT(DISTINCT ...)` guesswork.

**Melting the five cell-count columns into a long `cell_counts` table**
(instead of keeping `b_cell`, `cd8_t_cell`, ... as wide columns) is the
choice that matters most for scaling:

- Adding a new cell population (a new marker panel, a new assay) is an
  `INSERT` into `populations`, not an `ALTER TABLE` on a wide fact table —
  no migration, no downtime, no schema churn as the panel grows from 5
  populations to 50.
- Aggregate/analytics queries (relative frequency per sample, mean count
  per population per condition, population trends over
  `time_from_treatment_start`) become simple `GROUP BY population_id`
  queries instead of hand-written `UNION ALL`s over N wide columns.
- It keeps the row shape uniform regardless of how many populations exist,
  which matters once "hundreds of projects" means panels are no longer
  identical across projects — a wide table would force every project onto
  the same fixed column set or need per-project schema variants.

**Why this scales to hundreds of projects / thousands of samples /
mixed workloads:**

- All foreign keys are indexed (`idx_subjects_project`,
  `idx_samples_subject`, `idx_cell_counts_population`), so both
  "drill down" queries (project → subjects → samples → counts) and
  "roll up" analytics (aggregate counts by condition, treatment, or
  population across the whole database) hit indexes rather than full
  scans.
- Normalization keeps the write path (loading a new project's CSV) cheap
  and free of update anomalies, while the star-like shape of
  `cell_counts` (a fact table referencing `samples`/`populations`
  dimensions) is exactly the shape analytics/BI tools and SQL window
  functions are optimized for.
- SQLite comfortably handles millions of rows with this schema; if the
  project outgrows a single file (concurrent writers, multi-GB analytics
  workloads), the same schema ports directly to Postgres with no
  redesign.

## Code structure

Every analysis module follows the same shape: a function that takes a live
`sqlite3` connection and returns a `pandas.DataFrame`, plus a thin `main()`
CLI wrapper around it. `dashboard.py` imports those functions directly
instead of re-implementing the SQL/statistics, so the dashboard and the
standalone scripts can never disagree about a number — there is exactly
one place each query and each statistical test is written.

- **`load_data.py`** (Part 1) — reads `cell-count.csv` into a pandas
  `DataFrame`, validates it column-wise with vectorized boolean masks
  (missing columns, blank required fields, non-integer/negative numeric
  fields all reported with exact CSV line numbers before any database
  write happens), then drops/recreates the schema and loads it via
  `DataFrame.to_sql`. Wrapped in one transaction, so a bad CSV or an
  integrity conflict (duplicate sample id, a subject with conflicting
  metadata across rows) rolls back and leaves the previous database file
  untouched.
- **`analysis.py`** (Part 2) — `relative_frequencies(conn)`: one SQL query
  joining `cell_counts` → `samples` → `populations`, using a
  `SUM(...) OVER (PARTITION BY sample_id)` window function to get each
  sample's total in the same pass as its per-population counts, returning
  `sample, total_count, population, count, percentage`.
- **`response_analysis.py`** (Part 3) — compares miraclib responders vs.
  non-responders on baseline (`time_from_treatment_start = 0`) PBMC
  samples from melanoma subjects. Baseline-only is deliberate: each
  subject has exactly one day-0 sample, so filtering to it both (a) tests
  a *pre-treatment* biomarker rather than an effect of treatment/response
  itself, and (b) avoids pseudoreplication from treating three correlated
  draws per subject as independent observations. Runs a two-sided
  Mann-Whitney U test per population (percentages aren't guaranteed
  normal, so rank-based rather than a t-test) with Benjamini-Hochberg FDR
  correction across the five populations, and reports the rank-biserial
  effect size alongside p so a null result reads as "no effect" rather
  than "underpowered." `build_boxplot(df)` returns the Plotly figure
  (reused by the dashboard); `write_boxplot` saves it to
  `response_boxplots.html` for a standalone view.
- **`baseline_cohort.py`** (Part 4) — `baseline_cohort(conn)` pulls the
  same melanoma/miraclib/PBMC/baseline cohort as Part 3 (all responses,
  not just yes/no), and three small helpers roll it up by project, by
  response, and by sex.
- **`dashboard.py`** — a Streamlit app presenting Parts 2-4 against
  whatever `cell_counts.db` currently contains: the full relative-
  frequency table (filterable by sample), the response boxplots and
  stats table, and the Part 4 breakdowns. Streamlit was already a project
  dependency and turns a DataFrame-returning function into a browser
  table/chart in a couple of lines, so no separate frontend/API layer was
  needed for what is fundamentally "render these five functions' output."
- **`Makefile`** — `setup` / `pipeline` / `dashboard` wrap the
  install/build/serve steps above so the whole project runs without
  knowing (or needing) the underlying Python commands.

## Dashboard

```bash
make dashboard
```

Serves at **http://localhost:8501**. In GitHub Codespaces, the terminal
prints a "port forwarded" notification for `8501` (or open it from the
**Ports** tab) — click that link to open the dashboard in the browser;
Codespaces proxies it over HTTPS so no extra setup is needed.
