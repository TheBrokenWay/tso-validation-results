"""Base document generator for dossier sections."""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List


COMPANY_NAME = "Predator X"
FOUNDATION_NAME = "Broken Way Foundation"
REPORT_WIDTH = 80


def extract_compound_id(dossier: Dict[str, Any]) -> str:
    """Extract compound ID from dossier — checks all known locations.

    Raises ValueError if no compound ID found anywhere (data integrity issue).
    """
    sections = dossier.get("sections", {})
    candidate = dossier.get("candidate", {})
    locations = [
        dossier.get("compound_id"),
        # candidate.name only if it looks like a real ID (not "Genesis")
        candidate.get("name") if candidate.get("name", "").startswith("PRV_") else None,
        candidate.get("compound_id"),
        dossier.get("metadata", {}).get("compound_id"),
        sections.get("executive_summary", {}).get("compound_id"),
    ]
    for loc in locations:
        if loc and str(loc).strip() and loc != "UNKNOWN":
            return str(loc).strip()
    raise ValueError("Compound ID not found in dossier — data integrity issue")


def extract_disease_id(dossier: Dict[str, Any]) -> str:
    """Extract disease ID from dossier — checks all known locations.

    Raises ValueError if no disease ID found anywhere (data integrity issue).
    """
    sections = dossier.get("sections", {})
    candidate = dossier.get("candidate", {})
    ole = dossier.get("engines", {}).get("ole", {})
    ole_matched = ole.get("prv_diseases_matched", [])
    fin = dossier.get("finalization", {})
    disease_context = fin.get("disease_context", [])
    locations = [
        dossier.get("disease_id"),
        dossier.get("indication"),
        candidate.get("indication"),
        ole_matched[0] if ole_matched else None,
        disease_context[0] if isinstance(disease_context, list) and disease_context else None,
        sections.get("executive_summary", {}).get("disease_name"),
    ]
    for loc in locations:
        if loc and str(loc).strip() and loc not in ("UNKNOWN", "unknown"):
            return str(loc).strip()
    raise ValueError("Disease ID not found in dossier — data integrity issue")


class BaseDocumentGenerator:
    """Base class for all dossier document generators.

    Subclasses must set section_name and section_number, and implement generate().
    """

    section_name = "base"
    section_number = "00"

    def __init__(self, dossier_data: Dict[str, Any]):
        self.dossier = dossier_data
        self.sections_data = dossier_data.get("sections", {})

    def generate(self) -> str:
        """Generate formatted text report. Override in subclass."""
        raise NotImplementedError

    def _header(self, title: str) -> str:
        """Standard document header with compound/disease/tier metadata."""
        try:
            compound_id = extract_compound_id(self.dossier)
        except ValueError:
            compound_id = ""
        try:
            disease = extract_disease_id(self.dossier)
        except ValueError:
            disease = ""
        tier = self.dossier.get("tier", "") or self.dossier.get("finalization", {}).get("tier", "")
        lines = [
            "=" * REPORT_WIDTH,
            f"  {COMPANY_NAME} | {FOUNDATION_NAME}",
            f"  {title}",
            f"  Compound: {compound_id}" if compound_id else None,
            f"  Disease:  {disease}" if disease else None,
            f"  Tier:     {tier}" if tier else None,
            f"  Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
            "=" * REPORT_WIDTH,
            "",
        ]
        return "\n".join(l for l in lines if l is not None)

    def _section_header(self, title: str) -> str:
        """Section divider with title."""
        return f"\n{'─' * REPORT_WIDTH}\n  {title}\n{'─' * REPORT_WIDTH}\n"

    def _table(
        self,
        headers: List[str],
        rows: List[List[str]],
        col_widths: List[int] = None,
    ) -> str:
        """Render a formatted text table."""
        if not rows:
            return "  (no data)\n"
        if not col_widths:
            col_widths = [
                max(
                    len(str(h)),
                    max((len(str(r[i])) for r in rows), default=0),
                )
                + 2
                for i, h in enumerate(headers)
            ]
        lines = []
        header_line = "  ".join(
            str(h).ljust(w) for h, w in zip(headers, col_widths)
        )
        lines.append(f"  {header_line}")
        lines.append("  " + "  ".join("-" * w for w in col_widths))
        for row in rows:
            row_line = "  ".join(
                str(v).ljust(w) for v, w in zip(row, col_widths)
            )
            lines.append(f"  {row_line}")
        return "\n".join(lines) + "\n"

    def _metric_bar(self, value: float, target: float, width: int = 30) -> str:
        """Render a visual progress bar with pass/fail indicator."""
        pct = min(1.0, max(0.0, value / target)) if target > 0 else 0
        filled = int(pct * width)
        bar = "\u2588" * filled + "\u2591" * (width - filled)
        status = "PASS" if value <= target else "FAIL"
        return f"[{bar}] {value:.4f}/{target:.4f} {status}"

    def _kv(self, key: str, value: Any, indent: int = 2) -> str:
        """Key-value pair with indentation."""
        return f"{' ' * indent}{key}: {value}"

    def _footer(self) -> str:
        """Standard document footer with governance trace."""
        governance = self.sections_data.get("executive_summary", {}).get(
            "governance", {}
        )
        lines = [
            "",
            "\u2500" * REPORT_WIDTH,
            f"  CONFIDENTIAL -- {COMPANY_NAME} Proprietary",
            f"  Constitutional Compliance: Laws 1-51 verified",
            f"  Trace: {governance.get('trace_id', 'N/A')}",
            "=" * REPORT_WIDTH,
        ]
        return "\n".join(lines)
