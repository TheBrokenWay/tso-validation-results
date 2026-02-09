"""
Local patent index for FTO (Freedom-to-Operate). No live API calls.
Searchable by: patent numbers, title/abstract/claims text, InChIKey, SMILES, CAS, drug names, synonyms.
Refresh weekly via patent_data_refresh.py. Last refresh timestamp stored for audit.
"""
from __future__ import annotations

import json
import os
import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Optional
from datetime import datetime, timezone

# Default index path: PX_Data/patents/patent_index.db
_REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INDEX_PATH = _REPO_ROOT / "PX_Data" / "patents" / "patent_index.db"
METADATA_PATH = _REPO_ROOT / "PX_Data" / "patents" / "metadata.json"


def _get_index_path() -> Path:
    path = os.environ.get("PATENT_INDEX_PATH")
    return Path(path) if path else DEFAULT_INDEX_PATH


def _fts5_quote(text: Optional[str]) -> str:
    """Escape text for FTS5 MATCH so '=' and other operators don't cause syntax errors."""
    s = (text or "").strip()[:500]
    if not s:
        return ""
    return '"' + s.replace('"', '""') + '"'


def _ensure_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def init_index(index_path: Optional[Path] = None) -> Path:
    """Create schema if not exists. Returns path to index."""
    path = index_path or _get_index_path()
    _ensure_dir(path)
    conn = sqlite3.connect(str(path))
    try:
        conn.execute("""
            CREATE TABLE IF NOT EXISTS patents (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                patent_number TEXT,
                publication_number TEXT,
                application_number TEXT,
                title TEXT,
                abstract TEXT,
                claims TEXT,
                description TEXT,
                inventors TEXT,
                assignees TEXT,
                filing_date TEXT,
                grant_date TEXT,
                legal_status TEXT,
                family_id TEXT,
                source TEXT,
                created_at TEXT
            )
        """)
        conn.execute("CREATE INDEX IF NOT EXISTS idx_patents_grant_date ON patents(grant_date)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_patents_patent_number ON patents(patent_number)")
        conn.execute("""
            CREATE TABLE IF NOT EXISTS patent_chemicals (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                patent_id INTEGER NOT NULL,
                inchi_key TEXT,
                smiles TEXT,
                cas_number TEXT,
                chemical_name TEXT,
                drug_name TEXT,
                synonym TEXT,
                FOREIGN KEY (patent_id) REFERENCES patents(id)
            )
        """)
        conn.execute("CREATE INDEX IF NOT EXISTS idx_pc_patent_id ON patent_chemicals(patent_id)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_pc_inchi ON patent_chemicals(inchi_key)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_pc_smiles ON patent_chemicals(smiles)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_pc_drug_name ON patent_chemicals(drug_name)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_pc_cas ON patent_chemicals(cas_number)")
        # FTS5 for full-text search on title, abstract, claims (rowid = patents.id)
        conn.execute("""
            CREATE VIRTUAL TABLE IF NOT EXISTS patent_fts USING fts5(
                title,
                abstract,
                claims,
                content='patents',
                content_rowid='id'
            )
        """)
        # Triggers so FTS stays in sync with patents table
        conn.execute("""
            CREATE TRIGGER IF NOT EXISTS patents_ai AFTER INSERT ON patents BEGIN
                INSERT INTO patent_fts(rowid, title, abstract, claims) VALUES (new.id, new.title, new.abstract, new.claims);
            END
        """)
        conn.execute("""
            CREATE TRIGGER IF NOT EXISTS patents_ad AFTER DELETE ON patents BEGIN
                INSERT INTO patent_fts(patent_fts, rowid, title, abstract, claims) VALUES ('delete', old.id, old.title, old.abstract, old.claims);
            END
        """)
        conn.execute("""
            CREATE TRIGGER IF NOT EXISTS patents_au AFTER UPDATE ON patents BEGIN
                INSERT INTO patent_fts(patent_fts, rowid, title, abstract, claims) VALUES ('delete', old.id, old.title, old.abstract, old.claims);
                INSERT INTO patent_fts(rowid, title, abstract, claims) VALUES (new.id, new.title, new.abstract, new.claims);
            END
        """)
        conn.execute("""
            CREATE TABLE IF NOT EXISTS metadata (
                key TEXT PRIMARY KEY,
                value TEXT
            )
        """)
        conn.commit()
    finally:
        conn.close()
    return path


def get_last_refresh(index_path: Optional[Path] = None) -> Optional[datetime]:
    """Return last refresh UTC from metadata.json or DB."""
    path = index_path or _get_index_path()
    meta_file = path.parent / "metadata.json"
    if meta_file.exists():
        try:
            data = json.loads(meta_file.read_text(encoding="utf-8"))
            ts = data.get("last_refresh_utc")
            if ts:
                return datetime.fromisoformat(ts.replace("Z", "+00:00"))
        except Exception:
            pass
    db = path
    if db.exists():
        conn = sqlite3.connect(str(db))
        try:
            row = conn.execute("SELECT value FROM metadata WHERE key = ?", ("last_refresh_utc",)).fetchone()
            if row:
                return datetime.fromisoformat(row[0].replace("Z", "+00:00"))
        except Exception:
            pass
        finally:
            conn.close()
    return None


def set_last_refresh(utc: Optional[datetime] = None, index_path: Optional[Path] = None) -> None:
    """Write last refresh timestamp for audit."""
    path = index_path or _get_index_path()
    t = utc or datetime.now(timezone.utc)
    ts = t.isoformat()
    _ensure_dir(path)
    meta_file = path.parent / "metadata.json"
    data = {}
    if meta_file.exists():
        try:
            data = json.loads(meta_file.read_text(encoding="utf-8"))
        except Exception:
            pass
    data["last_refresh_utc"] = ts
    meta_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
    if path.exists():
        conn = sqlite3.connect(str(path))
        try:
            conn.execute(
                "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
                ("last_refresh_utc", ts),
            )
            conn.commit()
        finally:
            conn.close()


def search(
    smiles: Optional[str] = None,
    inchi_key: Optional[str] = None,
    drug_names: Optional[List[str]] = None,
    text_query: Optional[str] = None,
    cas_number: Optional[str] = None,
    limit: int = 100,
    index_path: Optional[Path] = None,
) -> List[Dict[str, Any]]:
    """
    Search local patent index. No external API calls.
    Returns list of patent dicts with patent_number, title, grant_date, assignees, etc.
    """
    path = index_path or _get_index_path()
    if not path.exists():
        return []
    conn = sqlite3.connect(str(path))
    conn.row_factory = sqlite3.Row
    try:
        seen_ids: set[int] = set()
        out: List[Dict[str, Any]] = []
        params: List[Any] = []
        clauses: List[str] = []

        # Chemical/drug exact match
        if smiles or inchi_key or cas_number or (drug_names and len(drug_names) > 0):
            sub = ["SELECT DISTINCT patent_id FROM patent_chemicals WHERE 1=0"]
            sub_params: List[Any] = []
            if smiles:
                sub.append(" OR smiles = ? OR smiles LIKE ?")
                sub_params.extend([smiles.strip()[:500], f"%{smiles.strip()[:200]}%"])
            if inchi_key:
                sub.append(" OR inchi_key = ? OR inchi_key LIKE ?")
                sub_params.extend([inchi_key.strip()[:100], f"%{inchi_key.strip()[:50]}%"])
            if cas_number:
                sub.append(" OR cas_number = ?")
                sub_params.append(cas_number.strip()[:50])
            if drug_names:
                placeholders = ",".join("?" * len(drug_names))
                sub.append(f" OR LOWER(drug_name) IN ({placeholders}) OR LOWER(synonym) IN ({placeholders})")
                dn = [d.strip().lower() for d in drug_names if d and len(d.strip()) >= 2][:20]
                sub_params.extend(dn)
                sub_params.extend(dn)
            if len(sub_params) > 0:
                subquery = "SELECT DISTINCT patent_id FROM patent_chemicals WHERE (1=1" + "".join(sub[1:]) + ")"
                clauses.append("p.id IN (" + subquery + ")")
                params.extend(sub_params)

        # Full-text search (FTS5 rowid = patents.id); quote literal to avoid syntax error on = etc.
        fts_query = _fts5_quote(text_query)
        if fts_query:
            fts_sql = "SELECT rowid FROM patent_fts WHERE patent_fts MATCH ? LIMIT ?"
            fts_ids = [row[0] for row in conn.execute(fts_sql, (fts_query, limit)).fetchall()]
            if fts_ids:
                placeholders = ",".join("?" * len(fts_ids))
                clauses.append(f"p.id IN ({placeholders})")
                params.extend(fts_ids)

        if not clauses:
            # No criteria: return nothing (don't dump whole DB)
            return []

        sql = "SELECT p.id, p.patent_number, p.publication_number, p.application_number, p.title, p.abstract, p.claims, p.inventors, p.assignees, p.filing_date, p.grant_date, p.legal_status, p.source FROM patents p WHERE " + " OR ".join(clauses) + " ORDER BY p.grant_date DESC LIMIT ?"
        params.append(limit)
        for row in conn.execute(sql, params).fetchall():
            pid = row["id"]
            if pid in seen_ids:
                continue
            seen_ids.add(pid)
            out.append({
                "patent_number": row["patent_number"] or row["publication_number"] or "",
                "title": (row["title"] or "")[:500],
                "abstract": (row["abstract"] or "")[:500],
                "grant_date": row["grant_date"] or "",
                "filing_date": row["filing_date"] or "",
                "assignees": row["assignees"] or "Unknown",
                "inventors": row["inventors"] or "",
                "legal_status": row["legal_status"] or "",
                "source": row["source"] or "",
            })
        return out
    finally:
        conn.close()


def index_exists(index_path: Optional[Path] = None) -> bool:
    """True if index file exists and has at least one patent."""
    path = index_path or _get_index_path()
    if not path.exists():
        return False
    conn = sqlite3.connect(str(path))
    try:
        n = conn.execute("SELECT COUNT(*) FROM patents").fetchone()[0]
        return n > 0
    except Exception:
        return False
    finally:
        conn.close()


def add_patent(
    patent_number: str,
    title: str = "",
    abstract: str = "",
    claims: str = "",
    description: str = "",
    inventors: str = "",
    assignees: str = "",
    filing_date: str = "",
    grant_date: str = "",
    legal_status: str = "",
    publication_number: str = "",
    application_number: str = "",
    family_id: str = "",
    source: str = "",
    chemicals: Optional[List[Dict[str, Any]]] = None,
    index_path: Optional[Path] = None,
) -> None:
    """Insert one patent and optional chemical rows. For use by refresh script."""
    path = index_path or _get_index_path()
    _ensure_dir(path)
    conn = sqlite3.connect(str(path))
    try:
        now = datetime.now(timezone.utc).isoformat()
        conn.execute(
            """INSERT INTO patents (
                patent_number, publication_number, application_number, title, abstract, claims, description,
                inventors, assignees, filing_date, grant_date, legal_status, family_id, source, created_at
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (
                patent_number or "",
                publication_number or "",
                application_number or "",
                title or "",
                abstract or "",
                claims or "",
                description or "",
                inventors or "",
                assignees or "",
                filing_date or "",
                grant_date or "",
                legal_status or "",
                family_id or "",
                source or "",
                now,
            ),
        )
        pid = conn.execute("SELECT last_insert_rowid()").fetchone()[0]
        # FTS is populated by AFTER INSERT trigger on patents
        for c in chemicals or []:
            conn.execute(
                """INSERT INTO patent_chemicals (patent_id, inchi_key, smiles, cas_number, chemical_name, drug_name, synonym)
                   VALUES (?,?,?,?,?,?,?)""",
                (
                    pid,
                    (c.get("inchi_key") or "")[:100],
                    (c.get("smiles") or "")[:500],
                    (c.get("cas_number") or "")[:50],
                    (c.get("chemical_name") or "")[:300],
                    (c.get("drug_name") or "")[:300],
                    (c.get("synonym") or "")[:300],
                ),
            )
        conn.commit()
    finally:
        conn.close()


__all__ = [
    "init_index",
    "get_last_refresh",
    "set_last_refresh",
    "search",
    "index_exists",
    "add_patent",
    "DEFAULT_INDEX_PATH",
    "METADATA_PATH",
]
