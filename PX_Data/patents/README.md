# Local Patent Dataset for FTO

Freedom-to-Operate checks use **only** this local searchable index. No live patent APIs are called.

## Data sources (public, legal)

- **USPTO Bulk Data** (downloaded by script): https://bulkdata.uspto.gov/ — patent grant full-text XML, weekly (Tuesdays). The refresh script downloads ZIPs from `data/patent/grant/redbook/fulltext/YYYY/`, parses XML, and indexes. Requires network access to bulkdata.uspto.gov. Use `--uspto-years 2024,2023` and `--uspto-max-zips N` (or env `PATENT_USPTO_YEARS`, `PATENT_USPTO_MAX_ZIPS`) to control how much to fetch.
- **USPTO ODP Bulk DataSets API (preferred, 2026)**: the ODP bulk directory is served via `api.uspto.gov` and **requires an ODP API key** (`X-API-KEY`). Set **`USPTO_ODP_API_KEY`** (or `USPTO_API_KEY`) in your **environment** (do not commit the key). Get a key at https://data.uspto.gov/apis/getting-started. The script discovers products/files and downloads via the Bulk DataSets API; use when `bulkdata.uspto.gov` does not resolve.
- **Lens Open Patent Data** (optional, with token): https://www.lens.org — bulk patent data as NDJSON.GZ. Set `LENS_ACCESS_TOKEN` (from Lens account → Subscriptions) to enable. Optionally set `PATENT_LENS_MAX_RECORDS` to cap records per run (full file is very large).
- **Google Patents Public Datasets** (no direct bulk HTTP): https://console.cloud.google.com/marketplace/product/google_patents_public_datasets — metadata, claims, citations via **BigQuery** only. To use: query BigQuery, export to GCS or JSONL, then load via `--seed` or a custom ingest.

## Refresh schedule

- **Every 7 days**: Run the refresh script to download new/updated patent files and rebuild the index.
- Incremental: re-download only new or updated files where the source provides manifests or timestamps.

## Storage layout

- `raw/` — downloaded bulk files (USPTO XML, etc.).
- `patent_index.db` — SQLite database with FTS5 full-text and chemical/drug indexes.
- `metadata.json` — last_refresh_utc, row counts, source URLs.

## FTO check procedure

- For each molecule/drug: convert to InChIKey and SMILES; search local index (exact + synonym + structure + drug name); return matching patents. No external API calls.
- Due diligence: weekly-refreshed local dataset satisfies reasonable FTO due diligence; refresh date and search results are logged for audit.

## Refresh schedule (7-day)

- Run the refresh script **every 7 days** (e.g. cron or Task Scheduler).
- On each run, re-download only new or updated patent files where the source provides manifests/timestamps; then re-index.
- Last refresh time is stored in `metadata.json` and in the DB table `metadata` for audit.

## Audit logging

- **Dataset refresh:** `metadata.json` in this directory holds `last_refresh_utc`. Each refresh run updates it.
- **FTO checks:** Every FTO check logs one line to `PX_LOGS/fto_audit.log` (repo root → PX_LOGS). Each line includes:
  - `timestamp` (UTC)
  - `compound_id`, `name`
  - `status` (CLEAR / BLOCKED / CLEAR_EXPIRED_PATENTS / LOCAL_INDEX_EMPTY / UNKNOWN)
  - `patent_hits`, `blocking` count
  - `dataset_refresh_utc` (index last refresh time)
- Retain `fto_audit.log` and `metadata.json` for due diligence and compliance.

## Commands

- **Refresh index (run weekly):**  
  `python -m PX_Executive.patent_data_refresh`  
  Or from repo root: `python PX_Executive/patent_data_refresh.py`  
  Options:
  - `--seed <path to JSONL>` — load from a seed file (one patent JSON per line).
  - `--uspto-years 2024,2023` — years to fetch from USPTO (default from `PATENT_USPTO_YEARS`).
  - `--uspto-max-zips N` — max ZIPs per year (default 4).
  - `--skip-uspto` — do not download USPTO; use seed and/or Lens only.
- **Check FTO:** Handled by `PX_Legal_Check`; it uses this index only (no live API calls).
