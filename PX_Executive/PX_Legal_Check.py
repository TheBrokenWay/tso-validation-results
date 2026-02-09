"""
╔══════════════════════════════════════════════════════════════════════════════╗
║ PX_Legal_Check.py                                                            ║
║ PREDATOR X :: PATENT FREEDOM TO OPERATE (LOCAL DATASET ONLY)                 ║
║ ARCHITECT: JAMES A. TILLAR | STATUS: ACTIVE | PROTOCOL: ZERO-COMPLIANT       ║
╚══════════════════════════════════════════════════════════════════════════════╝

PURPOSE:
    Freedom-to-operate (FTO) from a LOCAL patent dataset only. No live patent APIs are called.
    Data sources: USPTO Bulk Data, Google Patents Public Datasets, Lens Open Patent Data.
    Refresh: Run patent_data_refresh.py every 7 days. See PX_Data/patents/README.md.

FTO CHECK PROCEDURE:
    - For each molecule/drug: convert to InChIKey (from SMILES when available) and search
      the local index for exact/synonym/structure/drug-name matches.
    - Return all matching patents from the local dataset. Filter by grant date (20-year rule).
    - Due diligence: weekly-refreshed local dataset satisfies reasonable FTO due diligence.
    - Log dataset refresh date and FTO search results for audit.

FAIL-CLOSED:
    - If local index is missing or empty → BLOCK (status LOCAL_INDEX_EMPTY).
    - Unknown or no search terms → BLOCK when mode is REGULATORY.
"""
from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from datetime import datetime, timezone

# Optional RDKit for SMILES → InChIKey
try:
    from rdkit import Chem
    _HAS_RDKIT = True
except ImportError:
    _HAS_RDKIT = False


@dataclass
class PatentCheckResult:
    """Patent verification result."""
    freedom_to_operate: bool
    status: str  # "CLEAR", "BLOCKED", "CLEAR_EXPIRED_PATENTS", "LOCAL_INDEX_EMPTY", "UNKNOWN"
    patent_hits: int
    blocking_patents: List[Dict]
    search_terms: List[str]
    timestamp: str
    api_response_time_ms: float

    def to_dict(self) -> Dict:
        return {
            "freedom_to_operate": self.freedom_to_operate,
            "status": self.status,
            "patent_hits": self.patent_hits,
            "blocking_patents": self.blocking_patents,
            "search_terms": self.search_terms,
            "timestamp": self.timestamp,
            "api_response_time_ms": self.api_response_time_ms,
        }


def _smiles_to_inchi_key(smiles: str) -> Optional[str]:
    """Convert SMILES to InChIKey using RDKit if available."""
    if not _HAS_RDKIT or not smiles or len(smiles.strip()) < 2:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return None
        inchi = Chem.MolToInchi(mol)
        if not inchi:
            return None
        return Chem.MolToInchiKey(mol) or None
    except Exception:
        return None


class PX_Legal_Check:
    """
    FTO checker using ONLY the local patent index. No external API calls.
    Index is built and refreshed by patent_data_refresh.py (run every 7 days).
    """

    def __init__(self, api_key: Optional[str] = None, mode: str = "REGULATORY"):
        """
        Args:
            api_key: Ignored (kept for API compatibility).
            mode: "REGULATORY" (strict) or "RESEARCH" (permissive).
        """
        self.mode = mode
        self._index_path: Optional[Path] = None
        if os.environ.get("PATENT_INDEX_PATH"):
            self._index_path = Path(os.environ["PATENT_INDEX_PATH"])

    def _get_index_path(self) -> Path:
        if self._index_path is not None:
            return self._index_path
        from PX_Executive.Patent_Local_Index import DEFAULT_INDEX_PATH
        return DEFAULT_INDEX_PATH

    def _normalize_blocking(
        self,
        total_hits: int,
        blocking_patents: List[Dict],
        search_terms: List[str],
        timestamp: str,
        api_time_ms: float,
    ) -> PatentCheckResult:
        if blocking_patents:
            return PatentCheckResult(
                freedom_to_operate=False,
                status="BLOCKED",
                patent_hits=total_hits,
                blocking_patents=blocking_patents,
                search_terms=search_terms,
                timestamp=timestamp,
                api_response_time_ms=api_time_ms,
            )
        if total_hits > 0:
            return PatentCheckResult(
                freedom_to_operate=True,
                status="CLEAR_EXPIRED_PATENTS",
                patent_hits=total_hits,
                blocking_patents=[],
                search_terms=search_terms,
                timestamp=timestamp,
                api_response_time_ms=api_time_ms,
            )
        return PatentCheckResult(
            freedom_to_operate=True,
            status="CLEAR",
            patent_hits=0,
            blocking_patents=[],
            search_terms=search_terms,
            timestamp=timestamp,
            api_response_time_ms=api_time_ms,
        )

    def check_compound(
        self,
        smiles: str,
        iupac_name: Optional[str] = None,
        compound_id: Optional[str] = None,
        novel_invention: bool = False,
    ) -> PatentCheckResult:
        """
        Check compound against local patent index only. No API calls.
        Searches by SMILES, InChIKey (from SMILES), drug name, and text (title/abstract).

        novel_invention: If True, search ONLY by structure (SMILES/InChIKey). Do not search
        by name or text. Novel molecules we invent have no prior-art name; searching by
        label (e.g. "Genesis") would cause false hits from patents that mention that word.
        Result: no structure match → CLEAR (truly novel).
        """
        start_time = time.time()
        timestamp = datetime.now(timezone.utc).isoformat()
        search_terms: List[str] = []
        if not novel_invention and iupac_name:
            search_terms.append(iupac_name.strip())
        if smiles:
            search_terms.append(smiles[:50].strip())

        from PX_Executive.Patent_Local_Index import (
            index_exists,
            get_last_refresh,
            search as index_search,
            init_index,
        )

        index_path = self._get_index_path()
        init_index(index_path)

        if not index_exists(index_path):
            return PatentCheckResult(
                freedom_to_operate=False,
                status="LOCAL_INDEX_EMPTY",
                patent_hits=0,
                blocking_patents=[],
                search_terms=search_terms,
                timestamp=timestamp,
                api_response_time_ms=(time.time() - start_time) * 1000,
            )

        if not search_terms and not smiles:
            if self.mode == "REGULATORY":
                return PatentCheckResult(
                    freedom_to_operate=False,
                    status="UNKNOWN",
                    patent_hits=0,
                    blocking_patents=[],
                    search_terms=[],
                    timestamp=timestamp,
                    api_response_time_ms=(time.time() - start_time) * 1000,
                )
            return PatentCheckResult(
                freedom_to_operate=True,
                status="CLEAR",
                patent_hits=0,
                blocking_patents=[],
                search_terms=[],
                timestamp=timestamp,
                api_response_time_ms=(time.time() - start_time) * 1000,
            )

        inchi_key = _smiles_to_inchi_key(smiles) if smiles else None
        # Novel inventions: structure-only search (no name/text → avoids "Genesis" etc. false hits).
        if novel_invention:
            drug_names = None
            text_query = None
        else:
            drug_names = [t for t in search_terms if t and len(t) >= 2] if search_terms else []
            text_query = " ".join(search_terms[:3]) if search_terms else None

        hits = index_search(
            smiles=smiles[:500] if smiles else None,
            inchi_key=inchi_key,
            drug_names=drug_names[:10] if drug_names else None,
            text_query=text_query,
            limit=100,
            index_path=index_path,
        )
        api_time_ms = (time.time() - start_time) * 1000

        current_year = datetime.now().year
        blocking: List[Dict] = []
        for p in hits:
            date_str = (p.get("grant_date") or "")[:10]
            try:
                yr = int(date_str[:4]) if len(date_str) >= 4 else 0
                if yr and current_year - yr < 20:
                    blocking.append({
                        "patent_number": p.get("patent_number", ""),
                        "title": (p.get("title") or "")[:200],
                        "date": date_str,
                        "assignee": p.get("assignees", "Unknown"),
                        "years_remaining": 20 - (current_year - yr),
                    })
            except (ValueError, TypeError):
                pass

        result = self._normalize_blocking(len(hits), blocking, search_terms, timestamp, api_time_ms)

        # Audit: log refresh date and FTO result (sovereign log or file)
        try:
            last_refresh = get_last_refresh(index_path)
            refresh_str = last_refresh.isoformat() if last_refresh else "unknown"
            _log_fto_audit(
                compound_id=compound_id,
                smiles=smiles[:50] if smiles else None,
                name=iupac_name,
                status=result.status,
                patent_hits=result.patent_hits,
                blocking_count=len(result.blocking_patents),
                timestamp=timestamp,
                dataset_refresh_utc=refresh_str,
            )
        except Exception:
            pass

        return result

    def batch_check(self, compounds: List[Dict[str, str]]) -> List[PatentCheckResult]:
        """Check multiple compounds. No API rate limit (local only)."""
        return [
            self.check_compound(
                smiles=c.get("smiles", ""),
                iupac_name=c.get("iupac_name"),
                compound_id=c.get("compound_id"),
            )
            for c in compounds
        ]


def _log_fto_audit(
    compound_id: Optional[str],
    smiles: Optional[str],
    name: Optional[str],
    status: str,
    patent_hits: int,
    blocking_count: int,
    timestamp: str,
    dataset_refresh_utc: str,
) -> None:
    """Write FTO audit line for due diligence (dataset refresh date + result)."""
    try:
        from pathlib import Path
        repo = Path(__file__).resolve().parents[1]
        log_dir = repo / "PX_LOGS"
        log_dir.mkdir(parents=True, exist_ok=True)
        audit_file = log_dir / "fto_audit.log"
        line = (
            f"{timestamp}\t"
            f"compound_id={compound_id or ''}\t"
            f"name={name or ''}\t"
            f"status={status}\t"
            f"patent_hits={patent_hits}\t"
            f"blocking={blocking_count}\t"
            f"dataset_refresh_utc={dataset_refresh_utc}\n"
        )
        with open(audit_file, "a", encoding="utf-8") as f:
            f.write(line)
    except Exception:
        pass


# ═══════════════════════════════════════════════════════════════════════════
# PROTOCOL ZERO INTEGRATION
# ═══════════════════════════════════════════════════════════════════════════


def quick_check_for_protocol_zero(smiles: str, name: Optional[str] = None) -> Tuple[bool, Dict]:
    """Simplified patent check for Protocol Zero. Uses local index only."""
    checker = PX_Legal_Check(mode="RESEARCH")
    result = checker.check_compound(smiles=smiles, iupac_name=name)
    return result.freedom_to_operate, result.to_dict()


if __name__ == "__main__":
    print("═" * 80)
    print("PX_Legal_Check :: FTO (Local Patent Dataset Only)")
    print("═" * 80)

    from PX_Executive.Patent_Local_Index import index_exists, get_last_refresh, DEFAULT_INDEX_PATH
    if index_exists():
        refresh = get_last_refresh()
        print(f"\nLocal index: {DEFAULT_INDEX_PATH}")
        print(f"Last refresh: {refresh.isoformat() if refresh else 'unknown'}")
    else:
        print("\nLocal index empty or missing. Run: python PX_Executive/patent_data_refresh.py --seed <path to seed.jsonl>")
        print("See PX_Data/patents/README.md for data sources and 7-day refresh schedule.")

    test_smiles = "CC(=O)Oc1ccccc1C(=O)O"
    test_name = "Aspirin"
    checker = PX_Legal_Check(mode="RESEARCH")
    result = checker.check_compound(smiles=test_smiles, iupac_name=test_name)

    print(f"\nCompound: {test_name}")
    print(f"Freedom to Operate: {result.freedom_to_operate}")
    print(f"Status: {result.status}")
    print(f"Patent Hits: {result.patent_hits}")
    print(f"Blocking: {len(result.blocking_patents)}")
    print("═" * 80)
