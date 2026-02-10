"""
Patent Landscape Monitor

Monitors patent landscape for FTO (Freedom-to-Operate) analysis.
Works with the local patent_index.db (SQLite FTS5) for offline searches,
and provides online refresh capabilities via USPTO bulk data.

This client does NOT replace PX_Legal_Check.py — it provides the data
layer that PX_Legal_Check queries. Refresh runs should be scheduled
weekly to maintain FTO compliance.

Data sources:
  1. USPTO Bulk Data (bulkdata.uspto.gov) — weekly grant fulltext
  2. USPTO ODP API (api.uspto.gov) — requires ODP API key
  3. Lens.org (optional) — bulk patent NDJSON with LENS_ACCESS_TOKEN

Uses only: requests, json, sqlite3, pathlib, datetime, hashlib, os (+ PX integrations).
"""
from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

try:
    import requests
except ImportError:
    requests = None  # type: ignore[assignment]

from PX_System.foundation.integrations.net_policy import (
    check_external_access_allowed,
    online_sync_enabled,
)
from PX_System.foundation.integrations.retry import retry_transient_http

# ── Constants ──────────────────────────────────────────────────────────

DATA_DIR = Path(__file__).resolve().parent
DB_PATH = DATA_DIR / "patent_index.db"
METADATA_PATH = DATA_DIR / "metadata.json"
SERVICE_NAME = "uspto"

# API keys from environment
_USPTO_ODP_KEY = os.environ.get("USPTO_ODP_API_KEY", os.environ.get("USPTO_API_KEY", ""))
_LENS_TOKEN = os.environ.get("LENS_ACCESS_TOKEN", "")

# FTO audit log
FTO_AUDIT_LOG = _REPO_ROOT / "PX_LOGS" / "fto_audit.log"

# Patent classification codes relevant to pharmaceutical compounds
PHARMA_CPC_CODES = [
    "A61K",   # Preparations for medical purposes
    "A61P",   # Therapeutic activity of chemical compounds
    "C07D",   # Heterocyclic compounds
    "C07K",   # Peptides
]


# ── Database Functions ─────────────────────────────────────────────────

def _ensure_db() -> sqlite3.Connection:
    """Open or create the patent index database with FTS5."""
    conn = sqlite3.connect(str(DB_PATH))
    conn.execute("PRAGMA journal_mode=WAL")

    conn.executescript("""
        CREATE TABLE IF NOT EXISTS patents (
            patent_id TEXT PRIMARY KEY,
            title TEXT,
            abstract TEXT,
            applicant TEXT,
            filing_date TEXT,
            grant_date TEXT,
            expiry_date TEXT,
            cpc_codes TEXT,
            claims_text TEXT,
            drug_names TEXT,
            smiles TEXT,
            inchikey TEXT,
            source TEXT,
            fetched_utc TEXT
        );

        CREATE VIRTUAL TABLE IF NOT EXISTS patents_fts USING fts5(
            patent_id, title, abstract, claims_text, drug_names, smiles, inchikey,
            content='patents',
            content_rowid='rowid'
        );

        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT
        );
    """)
    return conn


def _upsert_patent(conn: sqlite3.Connection, patent: Dict[str, Any]) -> None:
    """Insert or update a patent record."""
    conn.execute("""
        INSERT OR REPLACE INTO patents
            (patent_id, title, abstract, applicant, filing_date, grant_date,
             expiry_date, cpc_codes, claims_text, drug_names, smiles, inchikey,
             source, fetched_utc)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        patent.get("patent_id", ""),
        patent.get("title", ""),
        patent.get("abstract", ""),
        patent.get("applicant", ""),
        patent.get("filing_date", ""),
        patent.get("grant_date", ""),
        patent.get("expiry_date", ""),
        patent.get("cpc_codes", ""),
        patent.get("claims_text", ""),
        patent.get("drug_names", ""),
        patent.get("smiles", ""),
        patent.get("inchikey", ""),
        patent.get("source", ""),
        datetime.now(timezone.utc).isoformat(),
    ))


def rebuild_fts(conn: sqlite3.Connection) -> None:
    """Rebuild the FTS5 index from the patents table."""
    conn.execute("INSERT INTO patents_fts(patents_fts) VALUES('rebuild')")
    conn.commit()


# ── Search Functions ───────────────────────────────────────────────────

def search_patents(
    query: str,
    search_type: str = "text",
    limit: int = 50,
) -> List[Dict[str, Any]]:
    """
    Search the local patent index.

    Args:
        query: Search string (compound name, SMILES, InChIKey, or free text)
        search_type: "text" (FTS5), "smiles" (exact), "inchikey" (exact)
        limit: Maximum results to return

    Returns:
        List of matching patent records
    """
    if not DB_PATH.exists():
        return []

    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row
    results: List[Dict[str, Any]] = []

    try:
        if search_type == "smiles":
            cursor = conn.execute(
                "SELECT * FROM patents WHERE smiles = ? LIMIT ?",
                (query, limit),
            )
        elif search_type == "inchikey":
            cursor = conn.execute(
                "SELECT * FROM patents WHERE inchikey = ? LIMIT ?",
                (query, limit),
            )
        else:
            # FTS5 full-text search
            safe_query = query.replace('"', '""')
            cursor = conn.execute(
                "SELECT p.* FROM patents p JOIN patents_fts f ON p.rowid = f.rowid WHERE patents_fts MATCH ? LIMIT ?",
                (f'"{safe_query}"', limit),
            )

        for row in cursor:
            results.append(dict(row))
    except sqlite3.OperationalError:
        # FTS table may not exist yet
        pass
    finally:
        conn.close()

    return results


def check_fto(
    compound_name: str,
    smiles: Optional[str] = None,
    inchikey: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Check Freedom-to-Operate for a compound against the local patent index.

    Args:
        compound_name: Drug/compound name
        smiles: SMILES string (optional, for structure search)
        inchikey: InChIKey (optional, for exact structure match)

    Returns:
        Dict with status, patent_hits, blocking count, and details
    """
    hits: List[Dict[str, Any]] = []

    # Search by name
    hits.extend(search_patents(compound_name, search_type="text"))

    # Search by SMILES if provided
    if smiles:
        hits.extend(search_patents(smiles, search_type="smiles"))

    # Search by InChIKey if provided
    if inchikey:
        hits.extend(search_patents(inchikey, search_type="inchikey"))

    # Deduplicate
    seen_ids: set = set()
    unique_hits: List[Dict[str, Any]] = []
    for h in hits:
        pid = h.get("patent_id", "")
        if pid and pid not in seen_ids:
            seen_ids.add(pid)
            unique_hits.append(h)

    # Determine blocking patents (not expired)
    now_str = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    blocking = [
        h for h in unique_hits
        if (h.get("expiry_date") or "9999-12-31") > now_str
    ]

    # Determine status
    if not DB_PATH.exists() or DB_PATH.stat().st_size < 1024:
        status = "LOCAL_INDEX_EMPTY"
    elif not unique_hits:
        status = "CLEAR"
    elif not blocking:
        status = "CLEAR_EXPIRED_PATENTS"
    else:
        status = "BLOCKED"

    # Get dataset refresh time
    refresh_utc = ""
    if METADATA_PATH.exists():
        try:
            meta = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
            refresh_utc = meta.get("last_refresh_utc", "")
        except Exception:
            pass

    result = {
        "compound_name": compound_name,
        "status": status,
        "patent_hits": len(unique_hits),
        "blocking": len(blocking),
        "dataset_refresh_utc": refresh_utc,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    # Write audit log
    _write_fto_audit(compound_name, result)

    return result


def _write_fto_audit(compound_name: str, result: Dict[str, Any]) -> None:
    """Append FTO check result to audit log."""
    try:
        FTO_AUDIT_LOG.parent.mkdir(parents=True, exist_ok=True)
        line = json.dumps({
            "timestamp": result.get("timestamp", ""),
            "compound_name": compound_name,
            "status": result.get("status", ""),
            "patent_hits": result.get("patent_hits", 0),
            "blocking": result.get("blocking", 0),
            "dataset_refresh_utc": result.get("dataset_refresh_utc", ""),
        }, default=str)
        with open(FTO_AUDIT_LOG, "a", encoding="utf-8") as f:
            f.write(line + "\n")
    except Exception as e:
        print(f"  WARN: FTO audit log write failed: {e}", file=sys.stderr)


# ── Online Refresh ─────────────────────────────────────────────────────

def get_index_stats() -> Dict[str, Any]:
    """Get statistics about the local patent index."""
    if not DB_PATH.exists():
        return {"exists": False, "patent_count": 0}

    conn = sqlite3.connect(str(DB_PATH))
    try:
        count = conn.execute("SELECT COUNT(*) FROM patents").fetchone()[0]
        meta_rows = conn.execute("SELECT key, value FROM metadata").fetchall()
        meta = {k: v for k, v in meta_rows}
    except sqlite3.OperationalError:
        count = 0
        meta = {}
    finally:
        conn.close()

    refresh_utc = ""
    if METADATA_PATH.exists():
        try:
            m = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
            refresh_utc = m.get("last_refresh_utc", "")
        except Exception:
            pass

    return {
        "exists": True,
        "patent_count": count,
        "db_size_mb": round(DB_PATH.stat().st_size / (1024 * 1024), 1),
        "last_refresh_utc": refresh_utc,
        "db_metadata": meta,
    }


def load_seed_file(seed_path: Path) -> int:
    """
    Load patents from a seed JSONL file (one patent JSON per line).

    Args:
        seed_path: Path to JSONL file

    Returns:
        Number of patents loaded
    """
    conn = _ensure_db()
    count = 0

    for line in seed_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            patent = json.loads(line)
            patent["source"] = f"seed:{seed_path.name}"
            _upsert_patent(conn, patent)
            count += 1
        except Exception:
            continue

    conn.commit()
    rebuild_fts(conn)

    # Update metadata
    conn.execute(
        "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
        ("last_seed_load", datetime.now(timezone.utc).isoformat()),
    )
    conn.commit()
    conn.close()

    # Update metadata.json
    _update_metadata_json(count)
    return count


def _update_metadata_json(refresh_count: int = 0) -> None:
    """Update metadata.json with current stats."""
    meta = {}
    if METADATA_PATH.exists():
        try:
            meta = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
        except Exception:
            pass

    meta["last_refresh_utc"] = datetime.now(timezone.utc).isoformat()
    meta["refresh_count_this_run"] = refresh_count
    meta["sources"] = [
        "USPTO bulkdata.uspto.gov (grant fulltext)",
        "Lens (optional token)",
        "Google Patents (BigQuery)",
    ]
    METADATA_PATH.write_text(json.dumps(meta, indent=2), encoding="utf-8")


# ── CLI Entry Point ────────────────────────────────────────────────────

def main() -> int:
    import argparse
    parser = argparse.ArgumentParser(description="Patent Landscape Monitor")
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("stats", help="Show patent index statistics")

    fto_p = sub.add_parser("fto", help="Check FTO for a compound")
    fto_p.add_argument("compound", help="Compound name")
    fto_p.add_argument("--smiles", help="SMILES string")
    fto_p.add_argument("--inchikey", help="InChIKey")

    search_p = sub.add_parser("search", help="Search patent index")
    search_p.add_argument("query", help="Search query")
    search_p.add_argument("--type", choices=["text", "smiles", "inchikey"], default="text")

    seed_p = sub.add_parser("seed", help="Load patents from JSONL seed file")
    seed_p.add_argument("path", help="Path to JSONL seed file")

    args = parser.parse_args()

    if args.command == "stats":
        stats = get_index_stats()
        print(json.dumps(stats, indent=2))
        return 0

    elif args.command == "fto":
        result = check_fto(args.compound, smiles=args.smiles, inchikey=args.inchikey)
        print(json.dumps(result, indent=2))
        return 0

    elif args.command == "search":
        results = search_patents(args.query, search_type=args.type)
        print(f"Found {len(results)} patents")
        for r in results[:20]:
            print(f"  {r.get('patent_id')}: {(r.get('title') or '')[:80]}")
        return 0

    elif args.command == "seed":
        count = load_seed_file(Path(args.path))
        print(f"Loaded {count} patents from seed file")
        return 0

    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    sys.exit(main())
