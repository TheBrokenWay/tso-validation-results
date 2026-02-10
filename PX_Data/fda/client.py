"""
FDA Guidance Feed Client

Monitors FDA drug approval and guidance data relevant to PRV-eligible diseases.
Uses two public APIs:
  1. openFDA Drug Approvals: https://api.fda.gov/drug/drugsfda.json
  2. FDA Guidance Documents: https://api.fda.gov/other/substance.json (and RSS feeds)

No API key required for basic usage (rate limited to 40 req/min without key).
Set FDA_API_KEY env var for higher rate limits.

Results are written to PX_Data/fda/ as JSONL + metadata.

Uses only: requests, json, csv, pathlib, datetime, hashlib, os (+ PX integrations).
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

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

OPENFDA_BASE = "https://api.fda.gov"
DATA_DIR = Path(__file__).resolve().parent
SERVICE_NAME = "fda"

# FDA API key (optional, for higher rate limits)
_API_KEY = os.environ.get("FDA_API_KEY", "")

# PRV-relevant search terms
PRV_DISEASE_TERMS = [
    "tropical disease", "neglected tropical disease",
    "nipah", "ebola", "marburg", "malaria", "tuberculosis",
    "chagas", "leishmaniasis", "dengue", "zika", "chikungunya",
    "rabies", "schistosomiasis", "cholera", "yellow fever",
    "priority review voucher", "rare pediatric disease",
    "orphan drug", "accelerated approval",
]

# Guidance document categories relevant to PRV pipeline
GUIDANCE_CATEGORIES = [
    "priority review", "tropical disease priority review voucher",
    "rare pediatric disease", "orphan drug designation",
    "accelerated approval", "breakthrough therapy",
]


# ── Core Functions ─────────────────────────────────────────────────────

def _build_params(**kwargs: Any) -> Dict[str, Any]:
    """Build request params, adding API key if available."""
    params = {k: v for k, v in kwargs.items() if v is not None}
    if _API_KEY:
        params["api_key"] = _API_KEY
    return params


@retry_transient_http(max_attempts=3, delay=2.0, backoff=2.0)
def _fetch_drug_approvals(search: str, limit: int = 100, skip: int = 0) -> Dict[str, Any]:
    """Fetch drug approval records from openFDA."""
    if requests is None:
        raise ImportError("requests library required for FDA API")

    params = _build_params(
        search=search,
        limit=limit,
        skip=skip,
    )
    resp = requests.get(f"{OPENFDA_BASE}/drug/drugsfda.json", params=params, timeout=30)
    resp.raise_for_status()
    return resp.json()


@retry_transient_http(max_attempts=3, delay=2.0, backoff=2.0)
def _fetch_drug_labels(search: str, limit: int = 100, skip: int = 0) -> Dict[str, Any]:
    """Fetch drug label records from openFDA."""
    if requests is None:
        raise ImportError("requests library required for FDA API")

    params = _build_params(
        search=search,
        limit=limit,
        skip=skip,
    )
    resp = requests.get(f"{OPENFDA_BASE}/drug/label.json", params=params, timeout=30)
    resp.raise_for_status()
    return resp.json()


def fetch_prv_approvals(max_results: int = 500) -> List[Dict[str, Any]]:
    """
    Fetch FDA drug approvals relevant to PRV-eligible diseases.

    Returns:
        List of approval records with key fields extracted
    """
    if not online_sync_enabled():
        print(f"  [fda] Network sync disabled", file=sys.stderr)
        return []
    if not check_external_access_allowed(SERVICE_NAME):
        print(f"  [fda] Access denied by net_policy", file=sys.stderr)
        return []

    all_approvals: List[Dict[str, Any]] = []

    for term in PRV_DISEASE_TERMS:
        try:
            search_q = f'openfda.brand_name:"{term}"+openfda.generic_name:"{term}"+products.active_ingredients.name:"{term}"'
            data = _fetch_drug_approvals(search=search_q, limit=min(100, max_results))
            results = data.get("results", [])

            for result in results:
                openfda = result.get("openfda", {})
                products = result.get("products", [])
                submissions = result.get("submissions", [])

                record = {
                    "application_number": result.get("application_number", ""),
                    "sponsor_name": result.get("sponsor_name", ""),
                    "brand_names": "|".join(openfda.get("brand_name", [])),
                    "generic_names": "|".join(openfda.get("generic_name", [])),
                    "substance_names": "|".join(openfda.get("substance_name", [])),
                    "product_count": len(products),
                    "submission_count": len(submissions),
                    "latest_submission_type": submissions[0].get("submission_type", "") if submissions else "",
                    "latest_submission_status": submissions[0].get("submission_status", "") if submissions else "",
                    "search_term": term,
                }

                # Extract active ingredients from products
                ingredients = []
                for p in products:
                    for ai in p.get("active_ingredients", []):
                        ingredients.append(ai.get("name", ""))
                record["active_ingredients"] = "|".join(set(ingredients))

                all_approvals.append(record)

            if results:
                print(f"  {term}: {len(results)} approvals")
            time.sleep(0.5)  # Rate limiting (40 req/min without key)

        except Exception as e:
            print(f"  {term}: FAILED ({e})", file=sys.stderr)

    return all_approvals


def fetch_prv_guidance() -> List[Dict[str, Any]]:
    """
    Fetch FDA guidance documents related to PRV regulatory pathway.

    Returns:
        List of guidance document summaries
    """
    if not online_sync_enabled():
        return []
    if not check_external_access_allowed(SERVICE_NAME):
        return []

    guidance_docs: List[Dict[str, Any]] = []

    for category in GUIDANCE_CATEGORIES:
        try:
            search_q = f'"{category}"'
            data = _fetch_drug_labels(search=search_q, limit=50)
            results = data.get("results", [])

            for result in results:
                openfda = result.get("openfda", {})
                doc = {
                    "brand_name": "|".join(openfda.get("brand_name", [])),
                    "generic_name": "|".join(openfda.get("generic_name", [])),
                    "indications_and_usage": (result.get("indications_and_usage") or [""])[0][:500],
                    "category": category,
                    "application_number": "|".join(openfda.get("application_number", [])),
                }
                guidance_docs.append(doc)

            if results:
                print(f"  guidance/{category}: {len(results)} docs")
            time.sleep(0.5)

        except Exception as e:
            print(f"  guidance/{category}: FAILED ({e})", file=sys.stderr)

    return guidance_docs


def refresh_all() -> Path:
    """
    Full refresh: fetch PRV approvals and guidance, write to JSONL + metadata.

    Returns:
        Path to the approvals output file
    """
    print(f"[fda] Refreshing FDA data for PRV pipeline...")

    # Fetch approvals
    approvals = fetch_prv_approvals()
    approvals_path = DATA_DIR / "fda_prv_approvals.jsonl"
    with open(approvals_path, "w", encoding="utf-8") as f:
        for record in approvals:
            f.write(json.dumps(record, default=str) + "\n")

    # Fetch guidance
    guidance = fetch_prv_guidance()
    guidance_path = DATA_DIR / "fda_prv_guidance.jsonl"
    with open(guidance_path, "w", encoding="utf-8") as f:
        for doc in guidance:
            f.write(json.dumps(doc, default=str) + "\n")

    # Write metadata
    meta = {
        "last_refresh_utc": datetime.now(timezone.utc).isoformat(),
        "approvals_count": len(approvals),
        "guidance_count": len(guidance),
        "source": "openFDA Drug Approvals + Labels API",
        "api_key_used": bool(_API_KEY),
        "checksum_approvals": hashlib.sha256(
            approvals_path.read_bytes() if approvals_path.exists() else b""
        ).hexdigest()[:16],
        "checksum_guidance": hashlib.sha256(
            guidance_path.read_bytes() if guidance_path.exists() else b""
        ).hexdigest()[:16],
    }
    meta_path = DATA_DIR / "metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"[fda] Done: {len(approvals)} approvals, {len(guidance)} guidance docs")
    return approvals_path


def search_approvals(compound_name: str) -> List[Dict[str, Any]]:
    """
    Search local JSONL for approvals matching a compound name.

    Args:
        compound_name: Drug/compound name to search for

    Returns:
        List of matching approval records
    """
    jsonl_path = DATA_DIR / "fda_prv_approvals.jsonl"
    if not jsonl_path.exists():
        return []

    matches: List[Dict[str, Any]] = []
    query_lower = compound_name.lower()

    for line in jsonl_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        searchable = " ".join([
            record.get("brand_names", ""),
            record.get("generic_names", ""),
            record.get("substance_names", ""),
            record.get("active_ingredients", ""),
        ]).lower()
        if query_lower in searchable:
            matches.append(record)

    return matches


# ── CLI Entry Point ────────────────────────────────────────────────────

def main() -> int:
    import argparse
    parser = argparse.ArgumentParser(description="FDA Guidance Feed client")
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("refresh", help="Refresh FDA approval and guidance data")

    search_p = sub.add_parser("search", help="Search local approvals by compound name")
    search_p.add_argument("compound", help="Compound name to search for")

    args = parser.parse_args()
    if args.command == "refresh":
        refresh_all()
        return 0
    elif args.command == "search":
        results = search_approvals(args.compound)
        print(f"Found {len(results)} matching approvals")
        for r in results[:20]:
            print(f"  {r.get('application_number')}: {r.get('brand_names', '')[:60]} ({r.get('generic_names', '')[:40]})")
        return 0
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    sys.exit(main())
