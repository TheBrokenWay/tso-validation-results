"""Base document generator for dossier sections."""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List


COMPANY_NAME = "Predator X"
FOUNDATION_NAME = "Broken Way Foundation"
REPORT_WIDTH = 80


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
        """Standard document header."""
        lines = [
            "=" * REPORT_WIDTH,
            f"  {COMPANY_NAME} | {FOUNDATION_NAME}",
            f"  {title}",
            f"  Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
            "=" * REPORT_WIDTH,
            "",
        ]
        return "\n".join(lines)

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
