"""
Patent local dataset refresh. Run every 7 days.
Downloads from USPTO Bulk Data (grant fulltext XML), optional Lens bulk (with token),
and documents Google Patents Public Datasets (BigQuery). Updates only new/changed
files where possible. Stores data in PX_Data/patents/ and builds the searchable index.
"""
from __future__ import annotations

import json
import os
import re
import sys
import tarfile
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime, timezone
from io import BytesIO
from typing import Any, Dict, List, Optional
import socket
import urllib.parse
import urllib.request

# Repo root
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

RAW_DIR = REPO_ROOT / "PX_Data" / "patents" / "raw"
INDEX_PATH = REPO_ROOT / "PX_Data" / "patents" / "patent_index.db"

# Data source URLs (public, legal)
#
# IMPORTANT (2026): `bulkdata.uspto.gov` is now frequently unreachable because the legacy BDSS
# DNS records have been retired in many environments. Use the ODP API (`api.uspto.gov`) instead
# for programmatic bulk dataset discovery + download.
USPTO_BULK_BASE = "https://bulkdata.uspto.gov"
USPTO_GRANT_FULLTEXT_BASE = f"{USPTO_BULK_BASE}/data/patent/grant/redbook/fulltext"
USPTO_ODP_API_BASE = "https://api.uspto.gov"
# Lens bulk: requires LENS_ACCESS_TOKEN (from Lens account Subscriptions tab)
LENS_BULK_RELEASE_URL = "https://api.lens.org/bulk/patent/release"
# Google Patents: BigQuery only — see README; no direct bulk HTTP.
GOOGLE_PATENTS_DOCS = "https://console.cloud.google.com/marketplace/product/google_patents_public_datasets"

# Default: download last 2 weeks of current year for USPTO (configurable)
USPTO_DEFAULT_YEARS = os.environ.get("PATENT_USPTO_YEARS", "2024")
USPTO_MAX_ZIPS = int(os.environ.get("PATENT_USPTO_MAX_ZIPS", "4"))


def ensure_dirs() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    INDEX_PATH.parent.mkdir(parents=True, exist_ok=True)


def _can_resolve(hostname: str) -> bool:
    """Return True if hostname resolves via system DNS."""
    try:
        socket.getaddrinfo(hostname, 443)
        return True
    except OSError:
        return False


def _get_uspto_odp_api_key() -> Optional[str]:
    """
    USPTO ODP API key (required for `api.uspto.gov` endpoints).
    Accepts several env var names for compatibility.
    """
    for k in ("USPTO_ODP_API_KEY", "USPTO_API_KEY", "ODP_API_KEY", "X_API_KEY"):
        v = os.environ.get(k, "").strip()
        if v:
            return v
    return None


def _http_get_odp_json(path: str, api_key: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    GET JSON from USPTO ODP API (api.uspto.gov) with X-API-KEY header.
    """
    url = USPTO_ODP_API_BASE.rstrip("/") + path
    if params:
        url += "?" + urllib.parse.urlencode({k: "" if v is None else v for k, v in params.items()}, doseq=True)
    req = urllib.request.Request(url, headers={"X-API-KEY": api_key, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = resp.read()
    return json.loads(data.decode("utf-8", errors="ignore") or "{}")


def _uspto_odp_search_products(api_key: str, q: str, limit: int = 25) -> List[Dict[str, Any]]:
    """
    Search bulk dataset products using ODP API.
    Endpoint: GET /api/v1/datasets/products/search?q=...&limit=...
    Response: { count, bulkDataProductBag: [ { productIdentifier, productTitleText, productFileBag... } ] }
    """
    payload = _http_get_odp_json("/api/v1/datasets/products/search", api_key, params={"q": q, "limit": limit})
    bag = payload.get("bulkDataProductBag") or payload.get("bulkDataProductBags") or []
    return bag if isinstance(bag, list) else []


def _uspto_odp_get_product(api_key: str, product_id: str, file_data_from_date: Optional[str] = None) -> Dict[str, Any]:
    """
    Fetch one product and its file listing.
    Endpoint: GET /api/v1/datasets/products/{productIdentifier}
    """
    params: Dict[str, Any] = {}
    if file_data_from_date:
        params["fileDataFromDate"] = file_data_from_date
    return _http_get_odp_json(f"/api/v1/datasets/products/{urllib.parse.quote(product_id)}", api_key, params=params)


def _uspto_odp_download_file(
    api_key: str,
    product_id: str,
    file_name: str,
    dest_path: Path,
    file_download_uri: Optional[str] = None,
) -> Path:
    """
    Download a bulk dataset file. Use file_download_uri if provided (path may include year);
    otherwise build URL as productId/fileName. urllib follows 302 redirects by default.
    """
    if file_download_uri and file_download_uri.startswith("http"):
        url = file_download_uri
    else:
        url = (
            USPTO_ODP_API_BASE.rstrip("/")
            + "/api/v1/datasets/products/files/"
            + urllib.parse.quote(product_id)
            + "/"
            + urllib.parse.quote(file_name)
        )
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"X-API-KEY": api_key, "User-Agent": "PredatorX-Patent-Refresh/1.0"})
    with urllib.request.urlopen(req, timeout=300) as resp, open(dest_path, "wb") as f:
        while True:
            chunk = resp.read(2**20)
            if not chunk:
                break
            f.write(chunk)
    return dest_path


def _http_get(url: str, timeout: int = 120, stream: bool = False):
    try:
        import requests
    except ImportError:
        from urllib.request import urlopen
        if stream:
            return urlopen(url, timeout=timeout)
        return urlopen(url, timeout=timeout).read()
    r = requests.get(url, timeout=timeout, stream=stream, headers={"User-Agent": "PredatorX-Patent-Refresh/1.0"})
    r.raise_for_status()
    if stream:
        return r
    return r.content


def _fetch_uspto_year_index(year: str) -> List[str]:
    """Return list of .zip filenames (e.g. ipg240102.zip) from USPTO fulltext year index page."""
    url = f"{USPTO_GRANT_FULLTEXT_BASE}/{year}/"
    try:
        data = _http_get(url, timeout=60)
        if hasattr(data, "read"):
            data = data.read()
        if isinstance(data, bytes):
            data = data.decode("utf-8", errors="ignore")
        # Prefer direct href="ipg240102.zip"; fallback: href=".../ipg240102.zip" or full URL
        direct = re.findall(r'href="(ipg\d{6}\.zip)"', data, re.I)
        if direct:
            return direct
        full = re.findall(r'href="[^"]*(ipg\d{6}\.zip)"', data, re.I)
        return list(dict.fromkeys(full))  # dedupe, preserve order
    except Exception:
        return []


def _tuesday_weeks_for_year(year: int, max_weeks: int = 52) -> List[str]:
    """Generate ipgYYMMDD.zip names for Tuesdays in year (first Tuesdays of year, so files likely exist)."""
    from datetime import date, timedelta
    out = []
    d = date(year, 1, 1)
    # Find first Tuesday
    while d.weekday() != 1:
        d += timedelta(days=1)
    while d.year == year and len(out) < max_weeks:
        out.append(f"ipg{d.strftime('%y%m%d')}.zip")
        d += timedelta(days=7)
    return out


def download_uspto_zip(year: str, zip_name: str, dest_dir: Path) -> Optional[Path]:
    """Download one USPTO weekly ZIP to dest_dir. Returns path to local ZIP or None."""
    url = f"{USPTO_GRANT_FULLTEXT_BASE}/{year}/{zip_name}"
    dest_dir.mkdir(parents=True, exist_ok=True)
    local = dest_dir / zip_name
    if local.exists():
        return local
    try:
        data = _http_get(url, timeout=300)
        if hasattr(data, "read"):
            data = data.read()
        local.write_bytes(data)
        return local
    except Exception as e:
        print(f"    Download failed {zip_name}: {e}", flush=True)
        return None


def _text(el: Optional[ET.Element]) -> str:
    if el is None:
        return ""
    return (el.text or "") + "".join((e.tail or "") for e in el).strip()


def _all_text(node: Optional[ET.Element]) -> str:
    if node is None:
        return ""
    return " ".join(node.itertext()).strip() if node.itertext else ""


def _local_tag(elem: ET.Element) -> str:
    """Return local tag name without namespace."""
    if elem.tag and "}" in elem.tag:
        return elem.tag.split("}", 1)[1]
    return elem.tag or ""


def _find_by_local_name(elem: ET.Element, local_name: str) -> Optional[ET.Element]:
    """Find first child or descendant with given local tag name."""
    for e in elem.iter():
        if _local_tag(e) == local_name:
            return e
    return None


def _findall_by_local_name(elem: ET.Element, local_name: str) -> List[ET.Element]:
    """Find all descendants with given local tag name."""
    return [e for e in elem.iter() if _local_tag(e) == local_name]


def parse_uspto_grant_xml_fragment(fragment: bytes) -> Optional[Dict[str, Any]]:
    """Parse one <us-patent-grant>...</us-patent-grant> XML fragment. Returns dict for normalize_patent."""
    try:
        root = ET.fromstring(fragment)
    except ET.ParseError:
        return None

    patent_number = ""
    publication_number = ""
    application_number = ""
    grant_date = ""
    filing_date = ""
    title = ""
    abstract = ""
    claims_text: List[str] = []
    description = ""
    inventors: List[str] = []
    assignees: List[str] = []

    bib = _find_by_local_name(root, "us-bibliographic-data-grant")
    if bib is None:
        bib = _find_by_local_name(root, "bibliographic-data")
    if bib is not None:
        pub_ref = _find_by_local_name(bib, "publication-reference")
        if pub_ref is not None:
            for doc_id in _findall_by_local_name(pub_ref, "document-id"):
                num = _find_by_local_name(doc_id, "doc-number")
                if num is not None:
                    publication_number = _text(num).strip()
                    patent_number = publication_number
                date_el = _find_by_local_name(doc_id, "date")
                if date_el is not None:
                    grant_date = _text(date_el).strip()[:10]
                break

        app_ref = _find_by_local_name(bib, "application-reference")
        if app_ref is not None:
            for doc_id in _findall_by_local_name(app_ref, "document-id"):
                anum = _find_by_local_name(doc_id, "doc-number")
                if anum is not None:
                    application_number = _text(anum).strip()
                date_el = _find_by_local_name(doc_id, "date")
                if date_el is not None:
                    filing_date = _text(date_el).strip()[:10]
                break

        for inv in _findall_by_local_name(bib, "inventor"):
            name = _find_by_local_name(inv, "name") or inv
            n = _all_text(name)
            if n:
                inventors.append(n)
        for a in _findall_by_local_name(bib, "assignee"):
            name = _find_by_local_name(a, "orgname")
            if name is None:
                name = _find_by_local_name(a, "name")
            if name is None:
                name = a
            n = _all_text(name)
            if n:
                assignees.append(n)

    tit = _find_by_local_name(root, "invention-title")
    if tit is not None:
        title = _all_text(tit)

    abs_node = _find_by_local_name(root, "abstract")
    if abs_node is not None:
        abstract = _all_text(abs_node)

    claims_node = _find_by_local_name(root, "claims")
    if claims_node is not None:
        for cl in _findall_by_local_name(claims_node, "claim"):
            claims_text.append(_all_text(cl))
    claims = " ".join(claims_text)[:20000]

    desc = _find_by_local_name(root, "description")
    if desc is not None:
        description = _all_text(desc)[:20000]

    return {
        "patent_number": patent_number or publication_number,
        "publication_number": publication_number,
        "application_number": application_number,
        "title": title,
        "abstract": abstract,
        "claims": claims,
        "description": description,
        "inventors": "; ".join(inventors)[:1000],
        "assignees": "; ".join(assignees)[:1000],
        "filing_date": filing_date,
        "grant_date": grant_date,
        "legal_status": "",
        "family_id": "",
        "source": "uspto",
        "chemicals": [],
    }


def extract_patents_from_uspto_zip(zip_path: Path) -> List[Dict[str, Any]]:
    """Extract and parse all patent grant XML documents from one USPTO weekly ZIP."""
    out: List[Dict[str, Any]] = []
    try:
        with zipfile.ZipFile(zip_path, "r") as z:
            for name in z.namelist():
                if not name.lower().endswith(".xml"):
                    continue
                with z.open(name) as f:
                    content = f.read()
                # Concatenated XML: multiple <us-patent-grant>...</us-patent-grant>
                start = content.find(b"<us-patent-grant")
                while start != -1:
                    end = content.find(b"</us-patent-grant>", start)
                    if end == -1:
                        break
                    end += len(b"</us-patent-grant>")
                    fragment = content[start:end]
                    rec = parse_uspto_grant_xml_fragment(fragment)
                    if rec and rec.get("patent_number"):
                        out.append(rec)
                    start = content.find(b"<us-patent-grant", end)
    except Exception:
        pass
    return out


def _parse_xml_chunk_for_grants(content: bytes) -> List[Dict[str, Any]]:
    """Split concatenated XML by <us-patent-grant> and parse each. Shared by zip/tar."""
    out: List[Dict[str, Any]] = []
    start = content.find(b"<us-patent-grant")
    while start != -1:
        end = content.find(b"</us-patent-grant>", start)
        if end == -1:
            break
        end += len(b"</us-patent-grant>")
        fragment = content[start:end]
        rec = parse_uspto_grant_xml_fragment(fragment)
        if rec and rec.get("patent_number"):
            out.append(rec)
        start = content.find(b"<us-patent-grant", end)
    return out


def extract_patents_from_uspto_tar(tar_path: Path) -> List[Dict[str, Any]]:
    """
    Extract and parse patent grant XML from one USPTO .tar (e.g. ODP PTGRDT weekly).
    PTGRDT tar contains per-patent .ZIP files; each ZIP typically contains one or more .xml with us-patent-grant.
    """
    out: List[Dict[str, Any]] = []
    try:
        with tarfile.open(tar_path, "r:*") as tar:
            for member in tar.getmembers():
                if not member.isfile():
                    continue
                name_lower = member.name.lower()
                try:
                    f = tar.extractfile(member)
                    if f is None:
                        continue
                    content = f.read()
                    if name_lower.endswith(".xml"):
                        out.extend(_parse_xml_chunk_for_grants(content))
                    elif name_lower.endswith(".zip"):
                        # Nested ZIP (one patent per ZIP in PTGRDT)
                        with zipfile.ZipFile(BytesIO(content), "r") as z:
                            for zname in z.namelist():
                                if not zname.lower().endswith(".xml"):
                                    continue
                                with z.open(zname) as zf:
                                    out.extend(_parse_xml_chunk_for_grants(zf.read()))
                except Exception:
                    continue
    except Exception:
        pass
    return out


def load_seed_jsonl(path: Path) -> List[Dict[str, Any]]:
    """Load patents from a JSONL file (one JSON object per line). For testing or manual ingest."""
    out = []
    if not path.exists():
        return out
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return out


def normalize_patent(rec: Dict[str, Any], source: str = "seed") -> Dict[str, Any]:
    """Normalize a patent record to our schema."""
    def s(key: str, default: str = "") -> str:
        v = rec.get(key)
        return (str(v).strip() if v is not None else default)[:2000]

    return {
        "patent_number": s("patent_number") or s("patent_id") or s("publication_number"),
        "publication_number": s("publication_number") or s("patent_number"),
        "application_number": s("application_number"),
        "title": s("title"),
        "abstract": s("abstract"),
        "claims": s("claims"),
        "description": s("description"),
        "inventors": s("inventors"),
        "assignees": s("assignees"),
        "filing_date": s("filing_date")[:10],
        "grant_date": s("grant_date")[:10],
        "legal_status": s("legal_status"),
        "family_id": s("family_id"),
        "source": rec.get("source") or source,
        "chemicals": rec.get("chemicals") or [],
    }


def download_and_index_uspto(
    years: str = USPTO_DEFAULT_YEARS,
    max_zips: int = USPTO_MAX_ZIPS,
    index_path: Optional[Path] = None,
) -> int:
    """
    Download USPTO patent grant fulltext ZIPs, parse XML, add to index.

    2026 reality check:
    - `bulkdata.uspto.gov` (legacy BDSS) often does not resolve in DNS any more.
    - Primary path is now the ODP Bulk DataSets API at `https://api.uspto.gov` which requires an API key.

    Control knobs:
    - Set `USPTO_ODP_API_KEY` (or `USPTO_API_KEY`) to enable ODP downloads.
    - Optional `PATENT_USPTO_PRODUCT_ID` to force a specific productIdentifier (shortName).
    - Optional `PATENT_USPTO_PRODUCT_QUERY` to search products (default: "Patent Grant Full Text").
    - If (and only if) bulkdata resolves, this function will also attempt legacy BDSS download
      by year/week (years + max_zips) as a fallback.

    Returns count of patents added.
    """
    from PX_Executive.Patent_Local_Index import init_index, add_patent
    path = index_path or INDEX_PATH
    init_index(path)
    added = 0

    # ──────────────────────────────────────────────────────────────────────────
    # Primary: ODP Bulk DataSets API (api.uspto.gov) — requires API key
    # ──────────────────────────────────────────────────────────────────────────
    api_key = _get_uspto_odp_api_key()
    product_id = os.environ.get("PATENT_USPTO_PRODUCT_ID", "").strip() or None
    product_query = os.environ.get("PATENT_USPTO_PRODUCT_QUERY", "Patent Grant Full Text").strip()

    if api_key:
        try:
            # Identify product
            product: Optional[Dict[str, Any]] = None
            if product_id:
                product = _uspto_odp_get_product(api_key, product_id)
            else:
                candidates = _uspto_odp_search_products(api_key, product_query, limit=25)
                # Prefer XML weekly products that actually provide downloadable zip files
                def score(p: Dict[str, Any]) -> int:
                    title = (p.get("productTitleText") or "").lower()
                    mimes = [str(x).upper() for x in (p.get("mimeTypeIdentifierArrayText") or [])]
                    freq = (p.get("productFrequencyText") or "").upper()
                    s = 0
                    if "XML" in mimes:
                        s += 10
                    if "grant" in title:
                        s += 5
                    if "full" in title and "text" in title:
                        s += 5
                    if freq == "WEEKLY":
                        s += 2
                    # Prefer items that have at least one downloadable fileName
                    file_bag = (p.get("productFileBag") or {})
                    file_data = file_bag.get("fileDataBag") if isinstance(file_bag, dict) else None
                    if isinstance(file_data, list) and any(isinstance(f, dict) and f.get("fileName") for f in file_data):
                        s += 3
                    return s

                candidates = sorted(candidates, key=score, reverse=True)
                product = candidates[0] if candidates else None

            if product and isinstance(product, dict):
                pid = product.get("productIdentifier") or product_id
                if pid:
                    # Get complete product details (includes fileDataBag)
                    # API returns { count, bulkDataProductBag: [ { productFileBag: { fileDataBag: [...] } } ] }
                    pfull = _uspto_odp_get_product(api_key, str(pid))
                    bag_list = (pfull.get("bulkDataProductBag") or []) if isinstance(pfull, dict) else []
                    first_product = bag_list[0] if isinstance(bag_list, list) and bag_list else {}
                    file_bag = first_product.get("productFileBag") if isinstance(first_product, dict) else None
                    files = (file_bag or {}).get("fileDataBag") if isinstance(file_bag, dict) else []
                    files = files if isinstance(files, list) else []

                    # Sort newest-first (release date or last modified), then take max_zips
                    def fkey(f: Dict[str, Any]) -> str:
                        return str(f.get("fileReleaseDate") or f.get("fileLastModifiedDateTime") or f.get("fileName") or "")

                    files = sorted([f for f in files if isinstance(f, dict) and f.get("fileName")], key=fkey, reverse=True)
                    files = files[:max_zips]

                    dest_dir = RAW_DIR / "uspto_odp" / str(pid)
                    print(f"USPTO ODP: product={pid} files={len(files)}", flush=True)
                    for f in files:
                        fname = str(f["fileName"])
                        # Preserve flat name for local path (URI may be productId/2026/file.tar)
                        local_name = fname.split("/")[-1] if "/" in fname else fname
                        local_file = dest_dir / local_name
                        if not local_file.exists():
                            print(f"  Downloading {fname} ...", flush=True)
                            _uspto_odp_download_file(
                                api_key, str(pid), fname, local_file,
                                file_download_uri=f.get("fileDownloadURI"),
                            )
                        else:
                            print(f"  Using cached {local_name}", flush=True)

                        # Parse .zip or .tar containing grant XML
                        suf = local_file.suffix.lower()
                        if suf == ".zip":
                            print(f"  Parsing {local_file.name} ...", flush=True)
                            records = extract_patents_from_uspto_zip(local_file)
                        elif suf == ".tar":
                            print(f"  Parsing {local_file.name} (tar) ...", flush=True)
                            records = extract_patents_from_uspto_tar(local_file)
                        else:
                            records = []
                        if records:
                            for rec in records:
                                norm = normalize_patent(rec, source="uspto_odp")
                                add_patent(
                                    index_path=path,
                                    patent_number=norm["patent_number"],
                                    title=norm["title"],
                                    abstract=norm["abstract"],
                                    claims=norm["claims"],
                                    description=norm["description"],
                                    inventors=norm["inventors"],
                                    assignees=norm["assignees"],
                                    filing_date=norm["filing_date"],
                                    grant_date=norm["grant_date"],
                                    legal_status=norm["legal_status"],
                                    publication_number=norm["publication_number"],
                                    application_number=norm["application_number"],
                                    family_id=norm["family_id"],
                                    source=norm["source"],
                                    chemicals=norm["chemicals"],
                                )
                                added += 1
                            print(f"    Indexed {len(records)} patents from {local_name}", flush=True)
                        else:
                            print(f"  Skipping (no grant XML parsed): {local_name}", flush=True)

                    return added
        except Exception as e:
            print(f"USPTO ODP download/index failed: {e}", flush=True)

    # ──────────────────────────────────────────────────────────────────────────
    # Legacy fallback: BDSS (bulkdata.uspto.gov) — only if it still resolves
    # ──────────────────────────────────────────────────────────────────────────
    if not _can_resolve("bulkdata.uspto.gov"):
        print(
            "USPTO BDSS host bulkdata.uspto.gov does not resolve in DNS here. "
            "Set USPTO_ODP_API_KEY to use the ODP Bulk DataSets API (api.uspto.gov).",
            flush=True,
        )
        return added

    dest_base = RAW_DIR / "uspto_bdss"
    year_list = [y.strip() for y in years.split(",") if y.strip()]
    for year in year_list:
        dest_dir = dest_base / year
        dest_dir.mkdir(parents=True, exist_ok=True)
        zip_names = _fetch_uspto_year_index(year)
        if not zip_names:
            try:
                y = int(year)
                zip_names = _tuesday_weeks_for_year(y, max_weeks=max_zips)
            except ValueError:
                continue
        for zip_name in zip_names[:max_zips]:
            local_zip = download_uspto_zip(year, zip_name, dest_dir)
            if local_zip is None:
                continue
            print(f"  Parsing {local_zip.name} ...", flush=True)
            records = extract_patents_from_uspto_zip(local_zip)
            for rec in records:
                norm = normalize_patent(rec, source="uspto_bdss")
                add_patent(
                    index_path=path,
                    patent_number=norm["patent_number"],
                    title=norm["title"],
                    abstract=norm["abstract"],
                    claims=norm["claims"],
                    description=norm["description"],
                    inventors=norm["inventors"],
                    assignees=norm["assignees"],
                    filing_date=norm["filing_date"],
                    grant_date=norm["grant_date"],
                    legal_status=norm["legal_status"],
                    publication_number=norm["publication_number"],
                    application_number=norm["application_number"],
                    family_id=norm["family_id"],
                    source=norm["source"],
                    chemicals=norm["chemicals"],
                )
                added += 1
            print(f"    Indexed {len(records)} patents from {zip_name}", flush=True)
    return added


def fetch_lens_release(token: str) -> Optional[Dict[str, Any]]:
    """Return latest Lens patent bulk release (requires LENS_ACCESS_TOKEN). Single object with downloadAccessKey, fileName, etc."""
    try:
        import requests
        r = requests.get(
            LENS_BULK_RELEASE_URL,
            headers={"Authorization": f"Bearer {token}", "User-Agent": "PredatorX-Patent-Refresh/1.0"},
            timeout=30,
        )
        r.raise_for_status()
        return r.json() if r.content else None
    except Exception:
        return None


def _lens_record_to_norm(rec: Dict[str, Any]) -> Dict[str, Any]:
    """Map one Lens patent JSON record to our normalized schema."""
    def s(key: str, default: str = "") -> str:
        v = rec.get(key)
        if v is None:
            return default
        if isinstance(v, dict):
            return (v.get("value") or v.get("text") or str(v))[:2000]
        return str(v).strip()[:2000]

    title = s("title")
    if not title and isinstance(rec.get("title"), list):
        titles = rec["title"]
        title = (titles[0].get("value") or titles[0].get("text") or "") if titles else ""

    pub_ref = rec.get("publication_reference") or rec.get("publicationReference") or {}
    pubs = pub_ref.get("publication") or pub_ref.get("publications") or []
    pub_number = ""
    pub_date = ""
    if isinstance(pubs, list) and pubs:
        p = pubs[0] if isinstance(pubs[0], dict) else {}
        pub_number = str(p.get("number") or p.get("doc_number") or rec.get("lens_id") or "")
        pub_date = (p.get("date") or "")[:10]
    if not pub_number:
        pub_number = str(rec.get("lens_id") or "")

    app_ref = rec.get("application_reference") or rec.get("applicationReference") or {}
    apps = app_ref.get("application") or app_ref.get("applications") or []
    app_number = ""
    filing_date = ""
    if isinstance(apps, list) and apps:
        a = apps[0] if isinstance(apps[0], dict) else {}
        app_number = str(a.get("number") or a.get("doc_number") or "")
        filing_date = (a.get("date") or "")[:10]

    abstract = s("abstract") or " ".join(
        (x.get("value") or x.get("text") or "") for x in (rec.get("abstract") or []) if isinstance(x, dict)
    )[:20000]
    claims = " ".join(
        (x.get("value") or x.get("text") or "") for x in (rec.get("claims") or []) if isinstance(x, dict)
    )[:20000]
    inventors = "; ".join(
        (x.get("name") or str(x)) for x in (rec.get("inventors") or []) if x
    )[:1000]
    assignees = "; ".join(
        (x.get("name") or str(x)) for x in (rec.get("assignees") or rec.get("applicants") or []) if x
    )[:1000]

    return {
        "patent_number": pub_number,
        "publication_number": pub_number,
        "application_number": app_number,
        "title": title,
        "abstract": abstract,
        "claims": claims,
        "description": "",
        "inventors": inventors,
        "assignees": assignees,
        "filing_date": filing_date,
        "grant_date": pub_date,
        "legal_status": "",
        "family_id": str(rec.get("family_id") or rec.get("familyId") or ""),
        "source": "lens",
        "chemicals": [],
    }


def download_and_index_lens(
    token: str,
    index_path: Optional[Path] = None,
    max_records: Optional[int] = None,
) -> int:
    """
    Fetch latest Lens patent bulk release, download .ndjson.gz, parse and add to index.
    max_records: cap for this run (None = process full file; full file is very large).
    Returns count of patents added.
    """
    import gzip

    from PX_Executive.Patent_Local_Index import init_index, add_patent

    path = index_path or INDEX_PATH
    init_index(path)
    release = fetch_lens_release(token)
    if not release or not release.get("downloadAccessKey"):
        print("Lens: no release or downloadAccessKey.", flush=True)
        return 0
    download_url = f"https://api.lens.org/bulk/download/{release['downloadAccessKey']}"
    dest_dir = RAW_DIR / "lens"
    dest_dir.mkdir(parents=True, exist_ok=True)
    local_file = dest_dir / (release.get("fileName") or "lens-patent-bulk.ndjson.gz")
    try:
        import requests
        r = requests.get(
            download_url,
            headers={"Authorization": f"Bearer {token}", "User-Agent": "PredatorX-Patent-Refresh/1.0"},
            stream=True,
            timeout=60,
        )
        r.raise_for_status()
        with open(local_file, "wb") as f:
            for chunk in r.iter_content(chunk_size=2**20):
                if chunk:
                    f.write(chunk)
    except Exception as e:
        print(f"Lens download failed: {e}", flush=True)
        return 0
    added = 0
    try:
        with gzip.open(local_file, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if max_records is not None and added >= max_records:
                    break
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                norm = _lens_record_to_norm(rec)
                if not norm.get("patent_number"):
                    continue
                add_patent(
                    index_path=path,
                    patent_number=norm["patent_number"],
                    title=norm["title"],
                    abstract=norm["abstract"],
                    claims=norm["claims"],
                    description=norm["description"],
                    inventors=norm["inventors"],
                    assignees=norm["assignees"],
                    filing_date=norm["filing_date"],
                    grant_date=norm["grant_date"],
                    legal_status=norm["legal_status"],
                    publication_number=norm["publication_number"],
                    application_number=norm["application_number"],
                    family_id=norm["family_id"],
                    source=norm["source"],
                    chemicals=norm["chemicals"],
                )
                added += 1
                if added % 10000 == 0:
                    print(f"  Lens: indexed {added} records ...", flush=True)
    except Exception as e:
        print(f"Lens parse error: {e}", flush=True)
    print(f"Lens: indexed {added} patents from {local_file.name}", flush=True)
    return added


def run_refresh(
    seed_path: Optional[Path] = None,
    incremental: bool = True,
    uspto_years: Optional[str] = None,
    uspto_max_zips: Optional[int] = None,
    skip_uspto: bool = False,
) -> None:
    """
    Refresh the local patent index.
    - If seed_path is set, load from that JSONL file (for testing or manual bulk).
    - If not skip_uspto, download USPTO grant fulltext ZIPs for uspto_years (default from env), parse and index.
    - Lens: set LENS_ACCESS_TOKEN to enable; then we fetch release list (download URLs need token + access key).
    - Google Patents: BigQuery only; see PX_Data/patents/README.md.
    """
    from PX_Executive.Patent_Local_Index import (
        init_index,
        add_patent,
        set_last_refresh,
        index_exists,
    )

    ensure_dirs()
    init_index(INDEX_PATH)
    count = 0

    if seed_path and seed_path.exists():
        records = load_seed_jsonl(seed_path)
        for rec in records:
            norm = normalize_patent(rec, source="seed")
            add_patent(
                index_path=INDEX_PATH,
                patent_number=norm["patent_number"],
                title=norm["title"],
                abstract=norm["abstract"],
                claims=norm["claims"],
                description=norm["description"],
                inventors=norm["inventors"],
                assignees=norm["assignees"],
                filing_date=norm["filing_date"],
                grant_date=norm["grant_date"],
                legal_status=norm["legal_status"],
                publication_number=norm["publication_number"],
                application_number=norm["application_number"],
                family_id=norm["family_id"],
                source=norm["source"],
                chemicals=norm["chemicals"],
            )
            count += 1
        print(f"Loaded {count} patents from seed file {seed_path}")

    if not skip_uspto:
        years = uspto_years or USPTO_DEFAULT_YEARS
        max_zips = uspto_max_zips if uspto_max_zips is not None else USPTO_MAX_ZIPS
        print(f"Downloading USPTO grant fulltext (years={years}, max_zips={max_zips}) ...")
        count += download_and_index_uspto(years=years, max_zips=max_zips, index_path=INDEX_PATH)

    # Lens: optional (requires LENS_ACCESS_TOKEN; full bulk is huge — set PATENT_LENS_MAX_RECORDS to cap)
    lens_token = os.environ.get("LENS_ACCESS_TOKEN", "").strip()
    lens_max = os.environ.get("PATENT_LENS_MAX_RECORDS", "")
    lens_max_int = int(lens_max) if lens_max.isdigit() else None
    if lens_token:
        try:
            print("Fetching and indexing Lens patent bulk (latest release) ...")
            count += download_and_index_lens(lens_token, index_path=INDEX_PATH, max_records=lens_max_int)
        except Exception as e:
            print(f"Lens bulk: {e}")

    set_last_refresh(datetime.now(timezone.utc), INDEX_PATH)
    meta = INDEX_PATH.parent / "metadata.json"
    data = {}
    if meta.exists():
        try:
            data = json.loads(meta.read_text(encoding="utf-8"))
        except Exception:
            pass
    data["last_refresh_utc"] = datetime.now(timezone.utc).isoformat()
    data["refresh_count_this_run"] = count
    data["sources"] = ["USPTO bulkdata.uspto.gov (grant fulltext)", "Lens (optional token)", "Google Patents (BigQuery)"]
    meta.write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(f"Last refresh set. Patents added this run: {count}")


def main() -> int:
    import argparse
    p = argparse.ArgumentParser(description="Refresh local patent dataset (run weekly)")
    p.add_argument("--seed", type=Path, default=None, help="Path to JSONL seed file (one patent JSON per line)")
    p.add_argument("--no-incremental", action="store_true", help="Full refresh (ignore incremental)")
    p.add_argument("--uspto-years", type=str, default=None, help=f"Comma-separated years (default: {USPTO_DEFAULT_YEARS})")
    p.add_argument("--uspto-max-zips", type=int, default=None, help=f"Max ZIPs per year to download (default: {USPTO_MAX_ZIPS})")
    p.add_argument("--skip-uspto", action="store_true", help="Do not download USPTO data (seed only)")
    args = p.parse_args()
    seed = args.seed
    if seed and not seed.is_absolute():
        seed = REPO_ROOT / seed
    run_refresh(
        seed_path=seed,
        incremental=not args.no_incremental,
        uspto_years=args.uspto_years,
        uspto_max_zips=args.uspto_max_zips,
        skip_uspto=args.skip_uspto,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
